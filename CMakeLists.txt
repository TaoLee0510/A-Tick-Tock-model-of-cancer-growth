cmake_minimum_required(VERSION 3.20)

project(ATCG LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS ON)

find_package(Freetype REQUIRED)
find_package(GSL REQUIRED)
find_package(PNG REQUIRED)
find_package(ZLIB REQUIRED)

if(APPLE AND CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    find_path(OPENMP_INCLUDE_DIR
        NAMES omp.h
        HINTS /opt/homebrew/opt/libomp/include /usr/local/opt/libomp/include
        REQUIRED
    )
    find_library(OPENMP_LIBRARY
        NAMES omp
        HINTS /opt/homebrew/opt/libomp/lib /usr/local/opt/libomp/lib
        REQUIRED
    )
    add_library(atcg_openmp INTERFACE)
    target_compile_options(atcg_openmp INTERFACE -Xclang -fopenmp)
    target_include_directories(atcg_openmp INTERFACE ${OPENMP_INCLUDE_DIR})
    target_link_libraries(atcg_openmp INTERFACE ${OPENMP_LIBRARY})
else()
    find_package(OpenMP REQUIRED)
    add_library(atcg_openmp INTERFACE)
    target_link_libraries(atcg_openmp INTERFACE OpenMP::OpenMP_CXX)
endif()

find_path(BLITZ_INCLUDE_DIR
    NAMES blitz/array.h
    HINTS /usr/local/include /opt/homebrew/include
    REQUIRED
)

find_path(PNGWRITER_INCLUDE_DIR
    NAMES pngwriter.h
    HINTS /usr/local/include /opt/homebrew/include
    REQUIRED
)

find_library(PNGWRITER_LIBRARY
    NAMES pngwriter
    HINTS /usr/local/lib /opt/homebrew/lib
    REQUIRED
)

add_executable(atcg
    ATCG/CountLines.cpp
    ATCG/cell_store.cpp
    ATCG/cell_trace.cpp
    ATCG/cell_motion.cpp
    ATCG/color_space.cpp
    ATCG/cell_type_transform.cpp
    ATCG/death_judgement.cpp
    ATCG/deltah_calculation.cpp
    ATCG/deltah_recalculation.cpp
    ATCG/density_calculation.cpp
    ATCG/density_growth_rate_calculation_1.cpp
    ATCG/division.cpp
    ATCG/free_living_division.cpp
    ATCG/free_living_division_single_thread.cpp
    ATCG/inner_count.cpp
    ATCG/inner_count_low_density.cpp
    ATCG/inner_initiation_array.cpp
    ATCG/int_grid.cpp
    ATCG/int_matrix.cpp
    ATCG/main.cpp
    ATCG/migrate_activation.cpp
    ATCG/migration.cpp
    ATCG/out_initiation_visualrange.cpp
    ATCG/outer_cell_count.cpp
    ATCG/outer_corr.cpp
    ATCG/outer_initiation_array.cpp
    ATCG/outer_initiation_array_low_density.cpp
    ATCG/read_files.cpp
    ATCG/random_migration.cpp
    ATCG/recovery_parameters.cpp
    ATCG/stateless_rng.cpp
    ATCG/stage_convert.cpp
    ATCG/visual_range.cpp
)

target_include_directories(atcg PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/ATCG
    ${BLITZ_INCLUDE_DIR}
    ${PNGWRITER_INCLUDE_DIR}
)

target_link_libraries(atcg PRIVATE
    Freetype::Freetype
    GSL::gsl
    GSL::gslcblas
    atcg_openmp
    PNG::PNG
    ZLIB::ZLIB
    ${PNGWRITER_LIBRARY}
)
