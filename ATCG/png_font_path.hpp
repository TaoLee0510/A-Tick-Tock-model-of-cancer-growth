#pragma once

#include <cstdlib>

inline const char* legacy_png_font_path() {
    const char* configured = std::getenv("ATCG_FONT_FILE");
    return configured && *configured ? configured : "Calisto MT.ttf";
}
