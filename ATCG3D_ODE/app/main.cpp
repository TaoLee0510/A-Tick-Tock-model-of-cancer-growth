#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>

#include "config/output_paths.hpp"
#include "ode/ode_model.hpp"

int main(int argc,char** argv) {
    try {
        std::filesystem::path path,root,report,checkpoint,resume;
        std::string model="ode";
        bool dry=false;
        for (int i=1;i<argc;++i) {
            const std::string key=argv[i];
            if (key=="--help") { std::cout<<"Usage: atcg3d_ode --config YAML [--model ode|periodic-pde] [--report JSON] [--output-root PATH]\n"
                <<"       [--checkpoint FILE] [--resume-checkpoint FILE] [--dry-run]\n"; return 0; }
            if (key=="--dry-run") { dry=true; continue; }
            if (++i>=argc) throw std::invalid_argument("missing argument value");
            if (key=="--config") path=argv[i]; else if(key=="--model") model=argv[i];
            else if(key=="--output-root") root=argv[i]; else if(key=="--report") report=argv[i];
            else if(key=="--checkpoint") checkpoint=argv[i]; else if(key=="--resume-checkpoint") resume=argv[i];
            else throw std::invalid_argument("unknown ODE option: "+key);
        }
        if (path.empty() || (model!="ode" && model!="periodic-pde")) throw std::invalid_argument("invalid ODE model/config");
        if(model!="ode" && (!checkpoint.empty() || !resume.empty())) throw std::invalid_argument("checkpoint options require --model ode");
        const auto config=atcg3d::ode::OdeConfig3D::load(path);
        if (dry) { std::cout<<"ODE configuration is valid\n"; return 0; }
        const auto directory=atcg3d::resolve_output_directory(config.output_directory,root);
        std::filesystem::create_directories(directory);
        std::ofstream metrics(directory/(resume.empty()?"metrics.csv":"metrics_resumed.csv"));
        if (!metrics) throw std::runtime_error("unable to write ODE metrics");
        metrics<<"time_hours,r_normal_small,r_normal_large,r_active_small,r_active_large,K_small,K_large,active_time_small,active_time_large,N,vessel_capacity,r_refractory_small,r_refractory_large,refractory_time_small,refractory_time_large\n"<<std::setprecision(17);
        atcg3d::ode::State final{};
        double time{};
        const auto row=[&](double t,const atcg3d::ode::State& s) { metrics<<t; for(double v:s) metrics<<','<<v; metrics<<'\n'; final=s; time=t; };
        if (model=="ode") {
            atcg3d::ode::OdeModel3D ode(config);
            if(!resume.empty()) ode.load_checkpoint(resume);
            row(ode.time_hours(),ode.state());
            while(ode.step()) row(ode.time_hours(),ode.state());
            if(!checkpoint.empty()) ode.save_checkpoint(checkpoint);
        } else {
            atcg3d::ode::PeriodicReactionPde3D pde(config); row(pde.time_hours(),pde.mean());
            while(pde.step()) row(pde.time_hours(),pde.mean());
        }
        if (report.empty()) report=directory/"summary.json";
        if (!report.parent_path().empty()) std::filesystem::create_directories(report.parent_path());
        std::ofstream out(report);
        if (!out) throw std::runtime_error("unable to write ODE summary");
        out<<std::setprecision(17)<<"{\"model\":\""<<model<<"\",\"time_hours\":"<<time<<",\"state\":[";
        for(std::size_t i=0;i<final.size();++i) { if(i) out<<','; out<<final[i]; }
        out<<"]}\n";
        std::cout<<"model="<<model<<" time_hours="<<time<<" mass_density="<<atcg3d::ode::Reaction3D(config.rules).total_mass(final)<<'\n';
        return 0;
    } catch(const std::exception& error) { std::cerr<<"ODE: "<<error.what()<<'\n'; return 1; }
}
