#pragma once

#include <cstddef>
#include <cstdint>
#include <atomic>
#include <chrono>
#include <condition_variable>
#include <deque>
#include <exception>
#include <filesystem>
#include <memory>
#include <mutex>
#include <optional>
#include <string>
#include <thread>
#include <vector>

#include "config/model_config.hpp"
#include "engine/simulation.hpp"
#include "io/run_manifest.hpp"

namespace atcg3d {

class OutputManager3D {
public:
    explicit OutputManager3D(const Model3DConfig& config);
    ~OutputManager3D() noexcept;

    OutputManager3D(const OutputManager3D&) = delete;
    OutputManager3D& operator=(const OutputManager3D&) = delete;

    void observe(const Simulation3D& simulation);
    void finalize(const Simulation3D& simulation);
    void checkpoint_now(const Simulation3D& simulation);

    const std::vector<SeriesEntry3D>& preview_entries() const noexcept { return preview_; }
    const std::vector<SeriesEntry3D>& full_entries() const noexcept { return full_; }
    const std::vector<SeriesEntry3D>& vessel_entries() const noexcept { return vessels_; }

private:
    struct FrozenPreviewState3D {
        std::vector<CellInit> cells;
        std::vector<Slot> cell_slots;
        std::vector<LesionId> lesion_ids;
        SimulationClock3D clock;
        SimulationStats3D stats;
        VasculatureState3D vasculature;
        std::size_t cell_slot_count{};
    };

    struct FrozenSimulationState3D {
        std::vector<CellInit> cells;
        std::vector<Slot> cell_slots;
        std::vector<Slot> free_slots;
        std::vector<LesionId> lesion_ids;
        CellUid next_uid{};
        SimulationClock3D clock;
        SimulationStats3D stats;
        std::vector<LineageEdge> lineage;
        VasculatureState3D vasculature;
        std::size_t cell_slot_count{};
        std::uint64_t state_checksum{};
    };

    struct FrozenCheckpointJournal3D {
        CheckpointCellJournal3D cells;
        CellUid next_uid{};
        SimulationClock3D clock;
        SimulationStats3D stats;
        std::vector<LineageEdge> lineage_tail;
        std::size_t lineage_prefix_count{};
        VasculatureState3D vasculature;
        std::uint64_t state_checksum{};
        std::filesystem::path path;
        std::filesystem::path parent_path;
        std::uint64_t parent_state_checksum{};
        double parent_time_hours{};
        std::uint64_t chain_length{};
    };

    struct OutputJob3D {
        std::optional<FrozenPreviewState3D> preview_state;
        std::optional<FrozenSimulationState3D> full_state;
        std::optional<FrozenCheckpointJournal3D> checkpoint_journal;
        bool preview{};
        bool live_preview{};
        bool full{};
        bool checkpoint_base{};
    };

    bool due(double now, double next, double interval) const noexcept;
    static double advance(double next, double now, double interval) noexcept;
    void write_preview(const Simulation3D& simulation);
    void write_preview(const FrozenPreviewState3D& simulation);
    void write_preview(const FrozenSimulationState3D& simulation);
    void write_live_preview(const FrozenPreviewState3D& simulation);
    void write_vessels(const Simulation3D& simulation, const std::string& frame);
    void write_vessels(const FrozenSimulationState3D& simulation,
                       const std::string& frame);
    void write_full(const Simulation3D& simulation);
    void write_full(const FrozenSimulationState3D& simulation);
    void write_checkpoint(const Simulation3D& simulation);
    void write_checkpoint(const FrozenSimulationState3D& simulation);
    void write_checkpoint(const FrozenCheckpointJournal3D& simulation);
    void append_lineage(const Simulation3D& simulation);
    void append_lineage(const std::vector<LineageEdge>& lineage);
    void validate_resume_state(const Simulation3D& simulation);
    void update_metadata();
    void write_metrics(const Simulation3D& simulation);
    void write_vascular_metrics(double time_hours,
                                const SimulationStats3D& stats,
                                const VasculatureState3D& vasculature);
    FrozenSimulationState3D freeze(const Simulation3D& simulation) const;
    FrozenPreviewState3D freeze_preview(const Simulation3D& simulation) const;
    FrozenCheckpointJournal3D freeze_checkpoint_journal(
        const Simulation3D& simulation,
        const std::filesystem::path& path);
    bool checkpoint_base_due(double now) const noexcept;
    void start_worker();
    static std::uint64_t estimate_job_bytes(const OutputJob3D& job) noexcept;
    void enqueue(OutputJob3D job);
    void worker_loop() noexcept;
    void finish_worker();
    void throw_worker_error();
    bool viewer_attached() const;

    Model3DConfig config_;
    std::filesystem::path run_directory_;
    std::vector<SeriesEntry3D> preview_;
    std::vector<SeriesEntry3D> full_;
    std::vector<SeriesEntry3D> vessels_;
    double next_preview_{};
    double next_full_{};
    double next_checkpoint_{};
    double last_preview_time_{-1.0};
    double last_full_time_{-1.0};
    double last_scheduled_preview_time_{-1.0};
    double last_scheduled_full_time_{-1.0};
    double last_scheduled_state_time_{-1.0};
    double resume_time_{-1.0};
    std::size_t lineage_written_{};
    std::vector<LineageEdge> existing_lineage_;
    std::filesystem::path checkpoint_parent_path_;
    std::filesystem::path last_checkpoint_path_;
    std::uint64_t checkpoint_parent_state_checksum_{};
    double checkpoint_parent_time_{-1.0};
    std::size_t checkpoint_lineage_scheduled_{};
    double checkpoint_base_time_{-1.0};
    std::uint64_t checkpoint_delta_chain_{};
    bool last_checkpoint_was_base_{true};
    bool resume_mode_{};
    bool resume_validated_{};
    bool finalized_{};
    std::mutex queue_mutex_;
    std::condition_variable queue_ready_;
    std::condition_variable queue_space_;
    std::deque<OutputJob3D> jobs_;
    std::uint64_t queued_snapshot_bytes_{};
    std::chrono::steady_clock::time_point last_live_preview_wall_{};
    std::thread worker_;
    std::exception_ptr worker_error_;
    std::atomic<bool> worker_failed_{false};
    bool worker_stopping_{};
};

}  // namespace atcg3d
