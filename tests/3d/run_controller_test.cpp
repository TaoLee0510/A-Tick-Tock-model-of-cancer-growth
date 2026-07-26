#include <cassert>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <string>
#include <thread>

#include "app/run_controller.hpp"

namespace {

std::string read_text(const std::filesystem::path& path) {
    std::ifstream stream(path, std::ios::binary);
    return {std::istreambuf_iterator<char>(stream),
            std::istreambuf_iterator<char>()};
}

void request_after(const std::filesystem::path& path,
                   std::chrono::milliseconds delay) {
    std::this_thread::sleep_for(delay);
    std::ofstream stream(path, std::ios::binary | std::ios::trunc);
    stream << "1\n";
}

}  // namespace

int main() {
    using namespace atcg3d;
    const std::filesystem::path directory =
        std::filesystem::temp_directory_path() /
        "atcg3d_run_controller_test";
    std::filesystem::remove_all(directory);

    Model3DConfig config;
    config.output_directory = directory;
    config.control_enabled = true;
    config.control_status_wall_interval_seconds = 0.05;
    config.control_poll_wall_interval_seconds = 0.01;
    config.end_time_hours = 10.0;
    config.initial_r_cells = 0;
    config.initial_K_cells = 0;
    config.angiogenesis.enabled = false;

    Simulation3D simulation(config);
    CellInit cell;
    cell.uid = 1;
    cell.clone_id = 1;
    cell.next_division_time = 10.0;
    cell.division_work_remaining = 10.0F;
    simulation.restore({cell}, 2, {}, {}, {});

    RunController3D controller(config);
    assert(controller.observe(simulation) ==
           RunControlAction3D::none);
    const auto control = directory / "control";
    assert(read_text(control / "status.json").find(
               "\"state\": \"running\"") != std::string::npos);

    {
        std::ofstream pause(control / "pause.request");
        pause << "1\n";
    }
    std::this_thread::sleep_for(std::chrono::milliseconds(15));
    std::thread resume(request_after, control / "resume.request",
                       std::chrono::milliseconds(80));
    const auto before = std::chrono::steady_clock::now();
    assert(controller.observe(simulation) ==
           RunControlAction3D::none);
    resume.join();
    assert(std::chrono::steady_clock::now() - before >=
           std::chrono::milliseconds(60));

    {
        std::ofstream pause(control / "pause.request");
        pause << "1\n";
    }
    std::this_thread::sleep_for(std::chrono::milliseconds(15));
    std::thread stop(request_after,
                     control / "checkpoint_stop.request",
                     std::chrono::milliseconds(40));
    assert(controller.observe(simulation) ==
           RunControlAction3D::checkpoint_and_stop);
    stop.join();
    controller.mark_stopped(simulation);
    assert(read_text(control / "status.json").find(
               "\"state\": \"stopped\"") != std::string::npos);

    std::filesystem::remove_all(directory);
    return 0;
}
