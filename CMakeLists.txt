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
    ATCG/cell_store.cpp
    ATCG/cell_trace.cpp
    ATCG/color_space.cpp
    ATCG/int_grid.cpp
    ATCG/int_matrix.cpp
    ATCG/main.cpp
    ATCG/recovery_parameters.cpp
    ATCG/stateless_rng.cpp
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
