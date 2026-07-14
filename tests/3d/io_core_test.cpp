#include <cassert>
#include <filesystem>
#include <fstream>
#include <string>

#include "config/model_config.hpp"
#include "core/cell_store.hpp"
#include "io/preview_sampler.hpp"
#include "io/run_manifest.hpp"

int main() {
    using namespace atcg3d;
    CellStore3D cells;
    for (std::uint64_t uid = 1; uid <= 100; ++uid) {
        CellInit cell;
        cell.uid = uid;
        cell.anchor = {static_cast<int>(uid), 0, 0};
        cells.create(cell);
    }
    const auto first = stable_preview_sample(cells, 10, 42);
    const auto second = stable_preview_sample(cells, 10, 42);
    const auto different = stable_preview_sample(cells, 10, 43);
    assert(first == second);
    assert(first.size() == 10);
    assert(first != different);
    assert(stable_preview_sample(cells, 200, 42).size() == 100);

    const auto directory = std::filesystem::temp_directory_path() / "atcg3d_manifest_test";
    std::filesystem::remove_all(directory);
    std::filesystem::create_directories(directory);
    const std::vector<SeriesEntry3D> preview{{"viz/preview/frame_00000000.vtkhdf", 0.0},
                                             {"viz/preview/frame_00000001.vtkhdf", 1.0}};
    const std::vector<SeriesEntry3D> full{{"viz/full/frame_00000000.vtkhdf", 0.0}};
    write_series_atomic(directory / "preview.vtkhdf.series", preview);
    Model3DConfig config;
    write_run_manifest_atomic(directory, config, preview, full, true);
    assert(std::filesystem::exists(directory / "preview.vtkhdf.series"));
    assert(std::filesystem::exists(directory / "run.json"));
    assert(!std::filesystem::exists(directory / "run.json.tmp"));
    std::ifstream stream(directory / "preview.vtkhdf.series");
    const std::string contents((std::istreambuf_iterator<char>(stream)), std::istreambuf_iterator<char>());
    assert(contents.find("file-series-version") != std::string::npos);
    assert(contents.find(".tmp") == std::string::npos);
    std::filesystem::remove_all(directory);
}
