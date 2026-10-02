#include "ode/ode_model.hpp"

#include <algorithm>
#include <bit>
#include <cmath>
#include <fstream>
#include <numeric>
#include <stdexcept>
#include <utility>

#include <yaml-cpp/yaml.h>
#include "common/density_growth_rule.hpp"

namespace atcg3d::ode {
namespace {
double growth_mean(const Model3DConfig& base, CellType type) {
    if (base.initial_growth_rate_model == "fixed") {
        return type == CellType::r ? base.initial_r_growth_rate : base.initial_K_growth_rate;
    }
    const auto& law = type == CellType::r ? base.initial_r_growth_truncated_normal : base.initial_K_growth_truncated_normal;
    return std::clamp(law.mean, law.minimum, law.maximum);
}
double beta_mean(const BetaRateConfig& law) { return law.scale * law.alpha / (law.alpha + law.beta); }
void check_state(const State& state) {
    for (const double value : state) if (!std::isfinite(value) || value < 0.0) throw std::invalid_argument("invalid ODE state");
    if (state[vessel_index] > 1.0 || state[10] > state[0] || state[11] > state[1] ||
        (state[2] == 0.0 && state[6] != 0.0) || (state[3] == 0.0 && state[7] != 0.0) ||
        (state[10] == 0.0 && state[12] != 0.0) || (state[11] == 0.0 && state[13] != 0.0)) {
        throw std::invalid_argument("inconsistent ODE compartment or clock");
    }
}
// Dormand-Prince embedded orders five and four. The error estimate controls
// each component; rejected steps never mutate the caller's state.
constexpr double a[7][7] = {
    {}, {1.0/5}, {3.0/40,9.0/40}, {44.0/45,-56.0/15,32.0/9},
    {19372.0/6561,-25360.0/2187,64448.0/6561,-212.0/729},
    {9017.0/3168,-355.0/33,46732.0/5247,49.0/176,-5103.0/18656},
    {35.0/384,0,500.0/1113,125.0/192,-2187.0/6784,11.0/84}
};
constexpr double b5[7] = {35.0/384,0,500.0/1113,125.0/192,-2187.0/6784,11.0/84,0};
constexpr double b4[7] = {5179.0/57600,0,7571.0/16695,393.0/640,-92097.0/339200,187.0/2100,1.0/40};
}  // namespace

OdeConfig3D OdeConfig3D::load(const std::filesystem::path& path) {
    const auto root = YAML::LoadFile(path.string());
    if (root["schema"]["name"].as<std::string>() != "atcg3d.ode_config" ||
        root["schema"]["version"].as<int>() != 1 ||
        root["model"].as<std::string>() != "well_mixed_structured_v1") throw std::invalid_argument("unsupported ODE schema/model");
    OdeConfig3D result;
    result.rules = structured_pde::StructuredPdeConfig3D::load(path.parent_path() / root["structured_config"].as<std::string>());
    if (const auto solver = root["solver"]) {
        if (solver["relative_tolerance"]) result.solver.relative_tolerance = solver["relative_tolerance"].as<double>();
        if (solver["absolute_tolerance"]) result.solver.absolute_tolerance = solver["absolute_tolerance"].as<double>();
        if (solver["maximum_step_hours"]) result.solver.maximum_step_hours = solver["maximum_step_hours"].as<double>();
    }
    const auto initial = root["initial"];
    for (const auto& [name, offset] : {std::pair{"r_normal",0}, {"r_active",2}, {"K",4}}) {
        const auto values = initial[name];
        if (!values.IsSequence() || values.size() != 2) throw std::invalid_argument("ODE populations require small/large pairs");
        for (int stage = 0; stage < 2; ++stage) result.initial[offset + stage] = values[stage].as<double>();
    }
    const auto clocks = initial["active_mean_hours"];
    if (!clocks.IsSequence() || clocks.size() != 2) throw std::invalid_argument("ODE requires two active mean clocks");
    for (int stage = 0; stage < 2; ++stage) result.initial[6 + stage] = clocks[stage].as<double>() * result.initial[2 + stage];
    result.initial[nutrient_index] = initial["nutrient"].as<double>();
    result.initial[vessel_index] = initial["vessel_capacity"].as<double>();
    result.output_directory = root["output"]["directory"].as<std::string>();
    result.validate();
    return result;
}
void OdeConfig3D::validate() const {
    rules.validate();
    if (rules.migration.activation_clock != "beta_mean_remaining_cycle_v1")
        throw std::invalid_argument("wrapper v1 does not carry activation duration distributions");
    if (rules.division_clock_model != "mean_rate_v1")
        throw std::invalid_argument("ODE v1 does not carry division work distributions");
    if (rules.schema_version < 7) throw std::invalid_argument("ODE requires the transported refractory contract v7+");
    for (const double value : {solver.relative_tolerance,solver.absolute_tolerance,solver.maximum_step_hours}) {
        if (!std::isfinite(value) || value <= 0.0) throw std::invalid_argument("invalid adaptive ODE tolerance/step");
    }
    check_state(initial);
    if (initial[nutrient_index] > rules.continuum.nutrient.vessel_value || output_directory.empty()) throw std::invalid_argument("invalid ODE initial nutrient/output");
    Reaction3D reaction(rules);
    if (reaction.occupied(initial) > rules.continuum.reaction.maximum_occupied_fraction) throw std::invalid_argument("ODE initial occupancy exceeds capacity");
}
Reaction3D::Reaction3D(structured_pde::StructuredPdeConfig3D rules) : rules_(std::move(rules)) {
    const auto& base = rules_.continuum.base;
    large_volume_ = base.thin_layer ? 4.0 : 8.0;
    window_measure_ = std::pow(static_cast<double>(base.growth_density_window_edge), base.thin_layer ? 2 : 3);
    r_inherent_ = growth_mean(base,CellType::r); K_inherent_ = growth_mean(base,CellType::K);
}
double Reaction3D::total_mass(const State& state) const noexcept {
    return std::accumulate(state.begin(),state.begin()+6,0.0);
}
double Reaction3D::occupied(const State& s) const noexcept {
    return s[0]+s[2]+s[4]+large_volume_*(s[1]+s[3]+s[5]);
}
State Reaction3D::derivative(const State& s, double rc, double kc) const {
    const auto& c = rules_.continuum;
    const auto& base = c.base;
    State d{};
    if (rc < 0.0) rc = (s[0]+s[1]+s[2]+s[3])*window_measure_;
    if (kc < 0.0) kc = (s[4]+s[5])*window_measure_;
    const double vacancy = std::clamp(1.0-occupied(s)/c.reaction.maximum_occupied_fraction,0.0,1.0);
    const double large_success = std::pow(vacancy,c.reaction.large_daughter_vacancy_exponent);
    const double small_success = 1.0-std::pow(1.0-vacancy,c.reaction.small_daughter_vacancy_exponent);
    const double density = total_mass(s);
    const auto growth = [&](CellType type, double inherent) {
        double rate = calculate_density_growth_rate_continuous(static_cast<int>(type),inherent,
            rc,kc,rc+kc,c.nutrient.common_density_limit,c.nutrient.common_density_limit,
            base.alpha,base.beta,c.nutrient.common_carrying_capacity,c.nutrient.common_carrying_capacity);
        return rate > 0.0 ? rate*s[8]/(c.nutrient.growth_half_saturation+s[8]) : rate;
    };
    const double gr = growth(CellType::r,r_inherent_), gk = growth(CellType::K,K_inherent_);
    const double lr = std::max(0.0,gr)/base.division_timing.base_cycle_hours;
    const double lk = std::max(0.0,gk)/base.division_timing.base_cycle_hours;
    const double dr = gr <= base.death_growth_rate_threshold ? 1.0/base.r_death_delay_hours : 0.0;
    const double dk = gk <= base.death_growth_rate_threshold ? 1.0/base.K_death_delay_hours : 0.0;
    const double conversion = base.r_to_K_conversion.enabled && density >= base.r_to_K_conversion.density_threshold
        ? base.r_to_K_conversion.probability_per_division : 0.0;
    const double small_loss = lr*(1.0-small_success)*c.reaction.failed_r_division_death_fraction;
    const double shape = lr*(1.0-large_success);
    for (int stage = 0; stage < 2; ++stage) {
        const double loss = dr + (stage == 0 ? small_loss : shape);
        d[stage] = -loss*s[stage]; d[2+stage] = -loss*s[2+stage];
        d[6+stage] = -s[2+stage]-loss*s[6+stage];
        d[10+stage] = -loss*s[10+stage];
        d[12+stage] = (s[12+stage] > 0.0 ? -s[10+stage] : 0.0)-loss*s[12+stage];
        d[4+stage] = -(dk+(stage == 1 ? lk*(1.0-large_success) : 0.0))*s[4+stage];
    }
    const double rsmall = s[0]+s[2], rlarge = s[1]+s[3];
    d[0] += (1.0-conversion)*small_success*lr*rsmall + (2.0-conversion)*shape*rlarge;
    d[1] += (1.0-conversion)*large_success*lr*rlarge;
    d[4] += conversion*(small_success*lr*rsmall+shape*rlarge) + small_success*lk*s[4] + 2.0*(1.0-large_success)*lk*s[5];
    d[5] += conversion*large_success*lr*rlarge + large_success*lk*s[5];
    const double rtotal = rsmall+rlarge, ktotal = s[4]+s[5];
    d[8] = -c.nutrient.decay_per_hour*s[8]
        -c.nutrient.r_consumption_rate_per_hour*rtotal*s[8]/(c.nutrient.r_consumption_half_saturation+s[8])
        -c.nutrient.K_consumption_rate_per_hour*ktotal*s[8]/(c.nutrient.K_consumption_half_saturation+s[8])
        +c.nutrient.vessel_exchange_per_hour*s[9]*(c.nutrient.vessel_value-s[8]);
    return d;
}
void Reaction3D::transitions(State& s) const {
    const double density = total_mass(s);
    for (int stage = 0; stage < 2; ++stage) {
        if (s[2+stage] > 0.0 && s[6+stage] <= 1.0e-12*std::max(1.0,s[2+stage])) {
            const double expired = s[2+stage];
            s[stage] += expired; s[10+stage] += expired;
            s[12+stage] += expired*rules_.migration.reactivation_cooldown_hours;
            s[2+stage] = s[6+stage] = 0.0;
        }
        if (s[12+stage] <= 1.0e-12 && density <= rules_.migration.reactivation_density_threshold) {
            s[10+stage] = s[12+stage] = 0.0;
        }
        if (density >= rules_.continuum.base.migration_activation_threshold) {
            const double mass = std::max(0.0,s[stage]-s[10+stage]);
            s[stage] -= mass; s[2+stage] += mass;
            s[6+stage] += mass*rules_.continuum.base.migration_activation_duration_mean_fraction*
                rules_.continuum.base.division_timing.base_cycle_hours/std::max(1.0e-12,r_inherent_);
        }
    }
}
void Reaction3D::advance(State& s, double duration, const SolverConfig& solver, double rc, double kc) const {
    if (!std::isfinite(duration) || duration < 0.0) throw std::invalid_argument("invalid ODE interval");
    check_state(s);
    double elapsed = 0.0, h = std::min(duration,solver.maximum_step_hours);
    for (std::size_t attempts = 0; elapsed < duration; ++attempts) {
        if (attempts > 1000000) throw std::runtime_error("adaptive ODE did not converge");
        transitions(s);
        h = std::min({h,solver.maximum_step_hours,duration-elapsed});
        for (int stage = 0; stage < 2; ++stage) {
            if (s[2+stage] > 0.0) h = std::min(h,s[6+stage]/s[2+stage]);
            if (s[10+stage] > 0.0 && s[12+stage] > 1.0e-12) h = std::min(h,s[12+stage]/s[10+stage]);
        }
        if (!(h > 0.0) || elapsed+h == elapsed) throw std::runtime_error("ODE step underflow");
        std::array<State,7> k{};
        for (int stage = 0; stage < 7; ++stage) {
            State work = s;
            for (std::size_t field = 0; field < s.size(); ++field)
                for (int previous = 0; previous < stage; ++previous) work[field] += h*a[stage][previous]*k[previous][field];
            k[stage] = derivative(work,rc,kc);
        }
        State candidate = s;
        double error = 0.0;
        bool positive = true;
        for (std::size_t field = 0; field < s.size(); ++field) {
            double estimate = 0.0;
            for (int stage = 0; stage < 7; ++stage) {
                candidate[field] += h*b5[stage]*k[stage][field];
                estimate += h*(b5[stage]-b4[stage])*k[stage][field];
            }
            error = std::max(error,std::abs(estimate)/(solver.absolute_tolerance+
                solver.relative_tolerance*std::max(std::abs(s[field]),std::abs(candidate[field]))));
            positive = positive && std::isfinite(candidate[field]) && candidate[field] >= -solver.absolute_tolerance;
        }
        if (error <= 1.0 && positive) {
            for (double& value : candidate) value = std::max(0.0,value);
            s = candidate; elapsed += h;
        }
        h *= positive ? std::clamp(error > 0.0 ? 0.9*std::pow(error,-0.2) : 5.0,0.1,5.0) : 0.5;
    }
    transitions(s);
    check_state(s);
}
OdeModel3D::OdeModel3D(OdeConfig3D config) : config_(std::move(config)),reaction_(config_.rules),state_(config_.initial),time_(config_.rules.continuum.start_time_hours) { config_.validate(); }
bool OdeModel3D::step() {
    const double end = config_.rules.continuum.end_time_hours;
    if (time_ >= end) return false;
    const double proposed = time_+config_.rules.continuum.time_step_hours;
    const double target = end-proposed <= 1.0e-12*std::max(1.0,end) ? end : proposed;
    reaction_.advance(state_,target-time_,config_.solver); time_ = target; return true;
}
void OdeModel3D::run() { while (step()) {} }
std::uint64_t OdeModel3D::state_checksum() const {
    std::uint64_t hash=config_.rules.dynamics_fingerprint();
    const auto mix=[&](double value) { hash^=std::bit_cast<std::uint64_t>(value); hash*=1099511628211ULL; };
    for(double value:{config_.solver.absolute_tolerance,config_.solver.relative_tolerance,config_.solver.maximum_step_hours,time_}) mix(value);
    for(double value:state_) mix(value);
    return hash;
}
void OdeModel3D::save_checkpoint(const std::filesystem::path& path) const {
    if(std::filesystem::exists(path)) throw std::runtime_error("refusing to overwrite ODE checkpoint");
    if(!path.parent_path().empty()) std::filesystem::create_directories(path.parent_path());
    std::ofstream out(path,std::ios::binary);
    const auto write=[&](const auto& value) { out.write(reinterpret_cast<const char*>(&value),sizeof(value)); };
    write(std::uint64_t{0x415443474f444531ULL}); write(config_.rules.dynamics_fingerprint());
    write(config_.solver.absolute_tolerance); write(config_.solver.relative_tolerance); write(config_.solver.maximum_step_hours);
    write(time_); write(state_); write(state_checksum());
    if(!out) throw std::runtime_error("unable to save ODE checkpoint");
}
void OdeModel3D::load_checkpoint(const std::filesystem::path& path) {
    std::ifstream in(path,std::ios::binary);
    const auto read=[&](auto& value) { in.read(reinterpret_cast<char*>(&value),sizeof(value)); if(!in) throw std::runtime_error("truncated ODE checkpoint"); };
    std::uint64_t magic{},fingerprint{},checksum{};
    double absolute{},relative{},maximum{},time{}; State state{};
    read(magic); read(fingerprint); read(absolute); read(relative); read(maximum); read(time); read(state); read(checksum);
    if(magic!=0x415443474f444531ULL || fingerprint!=config_.rules.dynamics_fingerprint() ||
        absolute!=config_.solver.absolute_tolerance || relative!=config_.solver.relative_tolerance || maximum!=config_.solver.maximum_step_hours ||
        !std::isfinite(time) || time<config_.rules.continuum.start_time_hours || time>config_.rules.continuum.end_time_hours ||
        in.peek()!=std::char_traits<char>::eof()) throw std::runtime_error("ODE checkpoint configuration mismatch");
    check_state(state);
    const auto previous_state=state_; const double previous_time=time_;
    state_=state; time_=time;
    if(state_checksum()!=checksum) { state_=previous_state; time_=previous_time; throw std::runtime_error("ODE checkpoint checksum mismatch"); }
}

PeriodicReactionPde3D::PeriodicReactionPde3D(OdeConfig3D config) : config_(std::move(config)),reaction_(config_.rules),time_(config_.rules.continuum.start_time_hours) {
    config_.validate(); const auto& shape = config_.rules.continuum.grid.shape;
    if(config_.rules.continuum.grid.spacing_voxels!=1.0) throw std::invalid_argument("periodic reaction reference requires unit spacing");
    fields_.assign(static_cast<std::size_t>(shape[0])*shape[1]*shape[2],config_.initial); work_ = fields_;
}
void PeriodicReactionPde3D::initialize(std::vector<State> fields) {
    if (fields.size() != fields_.size()) throw std::invalid_argument("periodic PDE shape mismatch");
    for (const auto& state : fields) check_state(state);
    fields_ = std::move(fields);
}
std::size_t PeriodicReactionPde3D::index(int x,int y,int z) const noexcept {
    const auto& shape = config_.rules.continuum.grid.shape;
    const auto wrap = [](int value,int count) { return (value%count+count)%count; };
    return (static_cast<std::size_t>(wrap(z,shape[2]))*shape[1]+wrap(y,shape[1]))*shape[0]+wrap(x,shape[0]);
}
State PeriodicReactionPde3D::mean() const {
    State result{};
    for (const auto& field : fields_) for (std::size_t i=0;i<result.size();++i) result[i] += field[i]/fields_.size();
    return result;
}
bool PeriodicReactionPde3D::step() {
    const auto& c = config_.rules.continuum; const auto& base = c.base;
    if (time_ >= c.end_time_hours) return false;
    const double proposed = time_+c.time_step_hours;
    const double target = c.end_time_hours-proposed <= 1.0e-12*std::max(1.0,c.end_time_hours) ? c.end_time_hours : proposed;
    const double dt = target-time_;
    const double factor = base.thin_layer ? 3.0/8 : 9.0/26;
    const double normal = factor*beta_mean(base.normal_r_migration_beta)*c.migration.diffusion_scale;
    const double Krate = factor*(base.initial_K_migration_rate_model == "fixed" ? base.initial_K_migration_rate : beta_mean(base.initial_K_migration_beta))*c.migration.diffusion_scale;
    const double active = normal*base.activated_r_normal_multiplier;
    const std::array<double,14> diffusivity{normal,normal,active,active,Krate,Krate,active,active,c.nutrient.diffusion_voxels2_per_hour,0,normal,normal,normal,normal};
    const int dimensions = base.thin_layer ? 2 : 3;
    const double maximum = *std::max_element(diffusivity.begin(),diffusivity.end());
    const int substeps = std::max(1,static_cast<int>(std::ceil(2.0*dimensions*dt*maximum/(0.45*c.grid.spacing_voxels*c.grid.spacing_voxels))));
    for (int substep = 0; substep < substeps; ++substep) {
        for (int z=0;z<c.grid.shape[2];++z) for (int y=0;y<c.grid.shape[1];++y) for (int x=0;x<c.grid.shape[0];++x) {
            const auto here=index(x,y,z);
            for (std::size_t f=0;f<14;++f) {
                double laplacian = fields_[index(x-1,y,z)][f]+fields_[index(x+1,y,z)][f]+fields_[index(x,y-1,z)][f]+fields_[index(x,y+1,z)][f]-4*fields_[here][f];
                if (!base.thin_layer) laplacian += fields_[index(x,y,z-1)][f]+fields_[index(x,y,z+1)][f]-2*fields_[here][f];
                work_[here][f] = fields_[here][f]+dt/substeps*diffusivity[f]*laplacian/(c.grid.spacing_voxels*c.grid.spacing_voxels);
            }
        }
        fields_.swap(work_);
    }
    const int edge = base.growth_density_window_edge, lower = (edge-1)/2, upper = edge-lower-1;
    for (int z=0;z<c.grid.shape[2];++z) for (int y=0;y<c.grid.shape[1];++y) for (int x=0;x<c.grid.shape[0];++x) {
        double rc=0,kc=0;
        for (int dz=base.thin_layer?0:-lower;dz<=(base.thin_layer?0:upper);++dz)
            for (int dy=-lower;dy<=upper;++dy) for (int dx=-lower;dx<=upper;++dx) {
                const auto& s=fields_[index(x+dx,y+dy,z+dz)]; rc+=s[0]+s[1]+s[2]+s[3]; kc+=s[4]+s[5];
            }
        const auto here=index(x,y,z);
        work_[here]=fields_[here]; reaction_.advance(work_[here],dt,config_.solver,rc,kc);
    }
    fields_.swap(work_); time_=target; return true;
}
}  // namespace atcg3d::ode
