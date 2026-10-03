#include "model/angiogenesis_field.hpp"
#include "model/shared_angiogenesis.hpp"

#include <algorithm>
#include <bit>
#include <cmath>
#include <iomanip>
#include <istream>
#include <numeric>
#include <ostream>
#include <sstream>
#include <stdexcept>

namespace atcg3d::continuum {
void AngiogenesisFieldConfig3D::validate() const {
    if(model!="disabled" && model!="vegf_tip_density_v1" && model!="shared_vegf_lattice_v2") throw std::invalid_argument("unknown angiogenesis field model");
    for(double v:{hypoxia_threshold,taf_production_per_cell_hour,taf_diffusion_voxels2_per_hour,taf_decay_per_hour,
        tip_diffusion_voxels2_per_hour,tip_chemotaxis,tip_branching_per_hour,tip_anastomosis_per_hour,
        seed_tips_per_hour,tip_speed_voxels_per_hour,vessel_radius_voxels,perfusion_exchange_per_hour,exclusion_fraction,maximum_tip_density}) {
        if(!std::isfinite(v) || v<0) throw std::invalid_argument("invalid angiogenesis field parameter");
    }
    if(hypoxia_threshold<=0 || hypoxia_threshold>1 || vessel_radius_voxels<=0 || maximum_tip_density<=0 ||
        exclusion_fraction<=0 || exclusion_fraction>1) throw std::invalid_argument("invalid vascular fraction/capacity");
}
std::string AngiogenesisFieldConfig3D::to_json() const {
    std::ostringstream out; out<<std::setprecision(17)<<"{\"model\":\""<<model<<"\"";
    const auto member=[&](const char* name,double value) { out<<",\""<<name<<"\":"<<value; };
    member("hypoxia_threshold",hypoxia_threshold); member("taf_production_per_cell_hour",taf_production_per_cell_hour);
    member("taf_diffusion_voxels2_per_hour",taf_diffusion_voxels2_per_hour); member("taf_decay_per_hour",taf_decay_per_hour);
    member("tip_diffusion_voxels2_per_hour",tip_diffusion_voxels2_per_hour); member("tip_chemotaxis",tip_chemotaxis);
    member("tip_branching_per_hour",tip_branching_per_hour); member("tip_anastomosis_per_hour",tip_anastomosis_per_hour);
    member("seed_tips_per_hour",seed_tips_per_hour); member("tip_speed_voxels_per_hour",tip_speed_voxels_per_hour);
    member("vessel_radius_voxels",vessel_radius_voxels); member("perfusion_exchange_per_hour",perfusion_exchange_per_hour);
    member("exclusion_fraction",exclusion_fraction); member("maximum_tip_density",maximum_tip_density);
    out<<'}'; return out.str();
}
std::uint64_t AngiogenesisFieldConfig3D::fingerprint() const {
    std::uint64_t h=1469598103934665603ULL; for(unsigned char ch:to_json()) { h^=ch; h*=1099511628211ULL; } return h;
}
AngiogenesisField3D::AngiogenesisField3D(AngiogenesisFieldConfig3D config,std::array<int,3> shape,double spacing,bool thin,int threads)
    :config_(std::move(config)),shape_(shape),spacing_(spacing),measure_(std::pow(spacing,thin?2:3)),
     cross_section_(thin?2*config_.vessel_radius_voxels:std::acos(-1.0)*config_.vessel_radius_voxels*config_.vessel_radius_voxels),thin_(thin) {
    config_.validate();
    if(spacing<=0 || !std::isfinite(spacing) || shape[0]<=0 || shape[1]<=0 || shape[2]<=0 || (thin && shape[2]!=1)) throw std::invalid_argument("invalid angiogenesis field grid");
    if (config_.model == "shared_vegf_lattice_v2") {
        shared_ = std::make_unique<SharedAngiogenesis3D>(config_, shape_, spacing_, thin_, threads);
        return;
    }
    const auto size=static_cast<std::size_t>(shape[0])*shape[1]*shape[2];
    for(auto* field:{&taf_,&tips_,&vessels_,&work_taf_,&work_tips_}) field->assign(size,0);
}
AngiogenesisField3D::~AngiogenesisField3D() = default;

void AngiogenesisField3D::use_individual_tips(std::uint64_t seed) {
    if (!shared_) throw std::logic_error("individual tips require the shared vascular model");
    shared_->use_individual_tips(seed);
}

const std::vector<double>& AngiogenesisField3D::taf() const noexcept {
    return shared_ ? shared_->taf() : taf_;
}

const std::vector<double>& AngiogenesisField3D::tips() const noexcept {
    return shared_ ? shared_->tips() : tips_;
}

const std::vector<double>& AngiogenesisField3D::vessels() const noexcept {
    return shared_ ? shared_->vessels() : vessels_;
}

const SharedVascularDiagnostics3D* AngiogenesisField3D::shared_diagnostics() const noexcept {
    return shared_ ? &shared_->diagnostics() : nullptr;
}

std::size_t AngiogenesisField3D::allocated_bytes() const noexcept {
    return shared_ ? shared_->allocated_bytes() : 5 * taf_.size() * sizeof(double);
}
std::size_t AngiogenesisField3D::index(int x,int y,int z) const noexcept { return (static_cast<std::size_t>(z)*shape_[1]+y)*shape_[0]+x; }
void AngiogenesisField3D::initialize(std::vector<double> vessels) {
    if (shared_) {
        shared_->initialize(std::move(vessels));
        return;
    }
    if(vessels.size()!=vessels_.size()) throw std::invalid_argument("vascular field shape mismatch");
    for(double value:vessels) if(!std::isfinite(value) || value<0 || value>1) throw std::invalid_argument("invalid vessel fraction");
    vessels_=std::move(vessels);
}
void AngiogenesisField3D::advance(double dt,const std::vector<double>& cells,const std::vector<double>& nutrient,double maximum,bool grow,
                                const VascularConsumerBounds3D* consumer_bounds) {
    if (shared_) {
        shared_->advance(dt, cells, nutrient, maximum, grow, consumer_bounds);
        return;
    }
    if(dt<0 || !std::isfinite(dt) || cells.size()!=taf_.size() || nutrient.size()!=taf_.size() || maximum<=0) throw std::invalid_argument("invalid vascular update");
    const int dims=thin_?2:3;
    std::vector<double> hypoxia(cells.size(),0),surface(cells.size(),0);
    double surface_sum=0,hypoxic_cells=0,total_cells=0;
    for(int z=0;z<shape_[2];++z) for(int y=0;y<shape_[1];++y) for(int x=0;x<shape_[0];++x) {
        const auto here=index(x,y,z);
        if(!std::isfinite(cells[here]) || cells[here]<0 || !std::isfinite(nutrient[here]) || nutrient[here]<0) throw std::invalid_argument("invalid vascular consumer/resource field");
        hypoxia[here]=std::clamp(1.0-nutrient[here]/(maximum*config_.hypoxia_threshold),0.0,1.0)*cells[here];
        total_cells+=cells[here]; if(nutrient[here]<maximum*config_.hypoxia_threshold) hypoxic_cells+=cells[here];
        const auto boundary=[&](int ax,int ay,int az) { return ax<0 || ax>=shape_[0] || ay<0 || ay>=shape_[1] || az<0 || az>=shape_[2] || cells[index(ax,ay,az)]<0.1*cells[here]; };
        if(cells[here]>0 && (boundary(x-1,y,z)||boundary(x+1,y,z)||boundary(x,y-1,z)||boundary(x,y+1,z)||
            (!thin_&&(boundary(x,y,z-1)||boundary(x,y,z+1))))) surface[here]=hypoxia[here];
        surface_sum+=surface[here];
    }
    if(surface_sum==0) { surface=hypoxia; surface_sum=std::accumulate(surface.begin(),surface.end(),0.0); }
    const double seed_rate=total_cells>0 ? config_.seed_tips_per_hour*hypoxic_cells/total_cells : 0;
    double elapsed=0;
    while(elapsed<dt) {
        double largest_rate=2.0*dims*std::max(config_.taf_diffusion_voxels2_per_hour,config_.tip_diffusion_voxels2_per_hour)/(spacing_*spacing_);
        for(int z=0;z<shape_[2];++z) for(int y=0;y<shape_[1];++y) for(int x=0;x<shape_[0];++x) {
            const auto here=index(x,y,z); double rate=2.0*dims*config_.tip_diffusion_voxels2_per_hour/(spacing_*spacing_);
            const auto face=[&](int ax,int ay,int az) { if(ax>=0&&ax<shape_[0]&&ay>=0&&ay<shape_[1]&&az>=0&&az<shape_[2]) rate+=config_.tip_chemotaxis*std::max(0.0,taf_[index(ax,ay,az)]-taf_[here])/(spacing_*spacing_); };
            face(x-1,y,z);face(x+1,y,z);face(x,y-1,z);face(x,y+1,z);if(!thin_){face(x,y,z-1);face(x,y,z+1);}
            largest_rate=std::max(largest_rate,rate);
        }
        const double h=std::min(dt-elapsed,largest_rate>0?0.45/largest_rate:dt-elapsed);
        if(!(h>0) || elapsed+h==elapsed) throw std::runtime_error("vascular CFL substep underflow");
        for(int z=0;z<shape_[2];++z) for(int y=0;y<shape_[1];++y) for(int x=0;x<shape_[0];++x) {
            const auto here=index(x,y,z); double lap=0,flux=0;
            const auto face=[&](int ax,int ay,int az) {
                if(ax<0||ax>=shape_[0]||ay<0||ay>=shape_[1]||az<0||az>=shape_[2]) return;
                const auto other=index(ax,ay,az); const double delta=taf_[other]-taf_[here]; lap+=delta;
                flux+=config_.tip_diffusion_voxels2_per_hour*(tips_[other]-tips_[here])/(spacing_*spacing_)
                    +config_.tip_chemotaxis*(std::max(0.0,-delta)*tips_[other]-std::max(0.0,delta)*tips_[here])/(spacing_*spacing_);
            };
            face(x-1,y,z);face(x+1,y,z);face(x,y-1,z);face(x,y+1,z);if(!thin_){face(x,y,z-1);face(x,y,z+1);}
            work_taf_[here]=std::max(0.0,(taf_[here]+h*config_.taf_diffusion_voxels2_per_hour*lap/(spacing_*spacing_)))*std::exp(-h*config_.taf_decay_per_hour)
                +h*config_.taf_production_per_cell_hour*hypoxia[here];
            if(!grow) { work_tips_[here]=0; continue; }
            const double moved=std::max(0.0,tips_[here]+h*flux);
            const double seeded=surface_sum>0 ? h*seed_rate*surface[here]/(surface_sum*measure_) : 0;
            const double branched=std::min(config_.maximum_tip_density,moved*std::exp(h*config_.tip_branching_per_hour*taf_[here]));
            work_tips_[here]=std::min(config_.maximum_tip_density,(branched+seeded)*std::exp(-h*config_.tip_anastomosis_per_hour*(vessels_[here]+moved)));
            vessels_[here]=1.0-(1.0-vessels_[here])*std::exp(-h*config_.tip_speed_voxels_per_hour*cross_section_*work_tips_[here]);
        }
        taf_.swap(work_taf_); tips_.swap(work_tips_); elapsed+=h;
    }
}
double AngiogenesisField3D::perfused_volume() const { return measure_*std::accumulate(vessels().begin(),vessels().end(),0.0); }
double AngiogenesisField3D::vessel_length() const { return shared_ ? shared_->diagnostics().centerline_growth : perfused_volume()/cross_section_; }
double AngiogenesisField3D::lesion_perfused_fraction(const std::vector<double>& cells) const {
    if(cells.size()!=vessels().size()) throw std::invalid_argument("lesion field shape mismatch");
    double lesion=0,perfused=0; for(std::size_t i=0;i<cells.size();++i) if(cells[i]>1.0e-12) { lesion+=measure_; perfused+=measure_*vessels()[i]; }
    return lesion>0?perfused/lesion:0;
}
std::uint64_t AngiogenesisField3D::checksum() const {
    if (shared_) return shared_->checksum();
    std::uint64_t hash=config_.fingerprint();for(const auto* field:{&taf_,&tips_,&vessels_}) for(double v:*field) { hash^=std::bit_cast<std::uint64_t>(v);hash*=1099511628211ULL; }return hash;
}
void AngiogenesisField3D::save(std::ostream& out) const {
    if (shared_) {
        shared_->save(out);
        return;
    }
    for(const auto* field:{&taf_,&tips_,&vessels_}) out.write(reinterpret_cast<const char*>(field->data()),field->size()*sizeof(double));
    if(!out) throw std::runtime_error("unable to save vascular fields");
}
void AngiogenesisField3D::load(std::istream& in) {
    if (shared_) {
        shared_->load(in);
        return;
    }
    for(auto* field:{&taf_,&tips_,&vessels_}) {
        in.read(reinterpret_cast<char*>(field->data()),field->size()*sizeof(double));
        for(double v:*field) if(!std::isfinite(v)||v<0||(field==&vessels_&&v>1)) throw std::runtime_error("invalid vascular checkpoint field");
    }
    if(!in) throw std::runtime_error("truncated vascular checkpoint");
}
}  // namespace atcg3d::continuum
