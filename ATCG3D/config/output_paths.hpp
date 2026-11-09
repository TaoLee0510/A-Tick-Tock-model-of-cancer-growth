#pragma once

#include <filesystem>

namespace atcg3d {

// An empty root retains YAML/CWD semantics. Explicit absolute output paths
// retain their meaning; only relative output directories are relocated.
inline std::filesystem::path resolve_output_directory(
    const std::filesystem::path& directory,
    const std::filesystem::path& root) {
    if (root.empty() || directory.is_absolute()) return directory;
    return (root / directory).lexically_normal();
}

}  // namespace atcg3d
