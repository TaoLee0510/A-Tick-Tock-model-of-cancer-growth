#include "app/run_controller.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <thread>

#if defined(__APPLE__) || defined(__linux__)
#include <unistd.h>
#endif

namespace atcg3d {
namespace {

template <class Writer>
void write_atomic(const std::filesystem::path& path, Writer writer) {
    std::filesystem::create_directories(path.parent_path());
    const std::filesystem::path temporary = path.string() + ".tmp";
    std::filesystem::remove(temporary);
    {
        std::ofstream stream(
            temporary, std::ios::binary | std::ios::trunc);
        if (!stream) {
            throw std::runtime_error(
                "unable to create control status " +
                temporary.string());
        }
        writer(stream);
        stream.flush();
        if (!stream) {
            throw std::runtime_error(
                "unable to finish control status " +
                temporary.string());
        }
    }
    std::filesystem::rename(temporary, path);
}

std::string json_escape(const std::string& value) {
    std::ostringstream out;
    for (const unsigned char character : value) {
        switch (character) {
            case '"': out << "\\\""; break;
            case '\\': out << "\\\\"; break;
            case '\n': out << "\\n"; break;
            case '\r': out << "\\r"; break;
            case '\t': out << "\\t"; break;
            default:
                if (character < 0x20U) {
                    out << "\\u" << std::hex << std::setw(4)
                        << std::setfill('0')
                        << static_cast<unsigned int>(character)
                        << std::dec;
                } else {
                    out << static_cast<char>(character);
                }
        }
    }
    return out.str();
}

long process_id() noexcept {
#if defined(__APPLE__) || defined(__linux__)
    return static_cast<long>(::getpid());
#else
    return 0;
#endif
}

}  // namespace

RunController3D::RunController3D(const Model3DConfig& config)
    : config_(config),
      control_directory_(config.output_directory / "control"),
      started_(std::chrono::steady_clock::now()),
      enabled_(config.control_enabled) {
    if (!enabled_) return;
    std::filesystem::create_directories(control_directory_);
    for (const char* name :
         {"pause.request", "resume.request",
          "checkpoint_stop.request", "terminate.request"}) {
        std::filesystem::remove(control_directory_ / name);
    }
    write_atomic(control_directory_ / "simulator.pid",
                 [](std::ostream& out) {
                     out << process_id() << '\n';
                 });
    publish(nullptr, "starting");
    // The first safe-point observation must immediately expose the live
    // simulation rather than leaving "starting" visible until the wall-clock
    // status interval expires.
    last_publish_ = {};
}

bool RunController3D::consume_request(const char* name) {
    const std::filesystem::path path = control_directory_ / name;
    if (!std::filesystem::exists(path)) return false;
    std::filesystem::remove(path);
    return true;
}

void RunController3D::publish(
    const Simulation3D* simulation,
    const char* state,
    const std::string& error) {
    if (!enabled_) return;
    const auto now = std::chrono::steady_clock::now();
    const double elapsed_seconds =
        std::chrono::duration<double>(now - started_).count();
    const double time_hours =
        simulation ? simulation->clock().time_hours : 0.0;
    const double progress =
        config_.end_time_hours > 0.0
            ? std::clamp(time_hours / config_.end_time_hours,
                         0.0, 1.0)
            : 1.0;
    const double simulated_hours_per_wall_hour =
        elapsed_seconds > 0.0
            ? time_hours * 3600.0 / elapsed_seconds
            : 0.0;
    const double eta_seconds =
        simulated_hours_per_wall_hour > 0.0
            ? (config_.end_time_hours - time_hours) /
                  simulated_hours_per_wall_hour * 3600.0
            : -1.0;
    write_atomic(control_directory_ / "status.json",
                 [&](std::ostream& out) {
        out << std::setprecision(17)
            << "{\n"
            << "  \"schema\": \"atcg3d.run-status\",\n"
            << "  \"schema_version\": 1,\n"
            << "  \"state\": \"" << state << "\",\n"
            << "  \"pid\": " << process_id() << ",\n"
            << "  \"time_hours\": " << time_hours << ",\n"
            << "  \"end_time_hours\": "
            << config_.end_time_hours << ",\n"
            << "  \"progress_fraction\": " << progress << ",\n"
            << "  \"elapsed_wall_seconds\": "
            << elapsed_seconds << ",\n"
            << "  \"simulated_hours_per_wall_hour\": "
            << simulated_hours_per_wall_hour << ",\n"
            << "  \"eta_wall_seconds\": " << eta_seconds << ",\n"
            << "  \"completed_events\": "
            << (simulation
                    ? simulation->clock().completed_events : 0)
            << ",\n"
            << "  \"alive_cells\": "
            << (simulation
                    ? simulation->cells().alive_count() : 0)
            << ",\n"
            << "  \"vessel_nodes\": "
            << (simulation
                    ? simulation->vessel_nodes().alive_count() : 0)
            << ",\n"
            << "  \"active_vessel_tips\": "
            << (simulation
                    ? simulation->active_vessel_tip_count() : 0)
            << ",\n"
            << "  \"threads\": " << config_.threads << ",\n"
            << "  \"output_directory\": \""
            << json_escape(config_.output_directory.string())
            << "\",\n"
            << "  \"error\": \"" << json_escape(error)
            << "\"\n"
            << "}\n";
    });
    last_publish_ = now;
}

void RunController3D::publish_if_due(
    const Simulation3D& simulation,
    const char* state) {
    const double elapsed =
        std::chrono::duration<double>(
            std::chrono::steady_clock::now() - last_publish_)
            .count();
    if (last_publish_.time_since_epoch().count() == 0 ||
        elapsed >=
            config_.control_status_wall_interval_seconds) {
        publish(&simulation, state);
    }
}

RunControlAction3D RunController3D::observe(
    const Simulation3D& simulation) {
    if (!enabled_ || terminal_) return RunControlAction3D::none;
    publish_if_due(simulation, "running");
    const auto now = std::chrono::steady_clock::now();
    if (last_poll_.time_since_epoch().count() != 0 &&
        std::chrono::duration<double>(now - last_poll_).count() <
            config_.control_poll_wall_interval_seconds) {
        return RunControlAction3D::none;
    }
    last_poll_ = now;
    if (consume_request("checkpoint_stop.request")) {
        publish(&simulation, "stopping");
        return RunControlAction3D::checkpoint_and_stop;
    }
    if (consume_request("terminate.request")) {
        publish(&simulation, "stopping");
        return RunControlAction3D::stop_without_checkpoint;
    }
    if (!consume_request("pause.request")) {
        return RunControlAction3D::none;
    }

    publish(&simulation, "paused");
    const auto delay = std::chrono::duration<double>(
        config_.control_poll_wall_interval_seconds);
    for (;;) {
        if (consume_request("checkpoint_stop.request")) {
            publish(&simulation, "stopping");
            return RunControlAction3D::checkpoint_and_stop;
        }
        if (consume_request("terminate.request")) {
            publish(&simulation, "stopping");
            return RunControlAction3D::stop_without_checkpoint;
        }
        if (consume_request("resume.request")) {
            last_poll_ = std::chrono::steady_clock::now();
            publish(&simulation, "running");
            return RunControlAction3D::none;
        }
        publish_if_due(simulation, "paused");
        std::this_thread::sleep_for(delay);
    }
}

void RunController3D::mark_completed(
    const Simulation3D& simulation) {
    terminal_ = true;
    publish(&simulation, "completed");
}

void RunController3D::mark_stopped(
    const Simulation3D& simulation) {
    terminal_ = true;
    publish(&simulation, "stopped");
}

void RunController3D::mark_failed(
    const std::string& message) noexcept {
    if (!enabled_ || terminal_) return;
    terminal_ = true;
    try {
        publish(nullptr, "failed", message);
    } catch (...) {
    }
}

}  // namespace atcg3d
