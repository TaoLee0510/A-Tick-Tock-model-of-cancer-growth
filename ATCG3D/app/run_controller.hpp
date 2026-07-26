#pragma once

#include <chrono>
#include <filesystem>
#include <string>

#include "config/model_config.hpp"
#include "engine/simulation.hpp"

namespace atcg3d {

enum class RunControlAction3D {
    none,
    checkpoint_and_stop,
    stop_without_checkpoint,
};

class RunController3D {
public:
    explicit RunController3D(const Model3DConfig& config);

    RunControlAction3D observe(const Simulation3D& simulation);
    void mark_completed(const Simulation3D& simulation);
    void mark_stopped(const Simulation3D& simulation);
    void mark_failed(const std::string& message) noexcept;

    const std::filesystem::path& control_directory() const noexcept {
        return control_directory_;
    }

private:
    bool consume_request(const char* name);
    void publish(const Simulation3D* simulation,
                 const char* state,
                 const std::string& error = {});
    void publish_if_due(const Simulation3D& simulation,
                        const char* state);

    Model3DConfig config_;
    std::filesystem::path control_directory_;
    std::chrono::steady_clock::time_point started_;
    std::chrono::steady_clock::time_point last_publish_{};
    std::chrono::steady_clock::time_point last_poll_{};
    bool enabled_{};
    bool terminal_{};
};

}  // namespace atcg3d
