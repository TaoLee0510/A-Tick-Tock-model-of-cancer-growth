cmake_minimum_required(VERSION 3.20)

project(ATCG LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS ON)

option(ATCG_BUILD_LEGACY_2D "Build the existing 2D ATCG executable" ON)
option(ATCG_BUILD_3D "Build the sparse event-driven 3D executable" ON)
option(ATCG3D_ENABLE_VTKHDF "Build VTK-HDF visualization output" OFF)
option(ATCG3D_ENABLE_HDF5_CHECKPOINT "Build HDF5 checkpoint support" OFF)
option(ATCG3D_ENABLE_LARGE_TESTS "Enable opt-in million/ten-million-cell tests" OFF)

if(ATCG_BUILD_LEGACY_2D OR ATCG_BUILD_3D)
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
endif()

if(ATCG_BUILD_LEGACY_2D)
    # These dependencies belong only to the historical 2D PNG-based target.
    find_package(Freetype REQUIRED)
    find_package(GSL REQUIRED)
    find_package(PNG REQUIRED)
    find_package(ZLIB REQUIRED)

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
        ATCG/CellDivision.cpp
        ATCG/CellDivisionSingleCell.cpp
        ATCG/CellMigration.cpp
        ATCG/CellMigrationDivision.cpp
        ATCG/CellMigrationDivisionSingleCell.cpp
        ATCG/CountLines.cpp
        ATCG/SaveAllPNG.cpp
        ATCG/SaveCellArray.cpp
        ATCG/SaveCellArraySingleCell.cpp
        ATCG/SaveCellTraceArray.cpp
        ATCG/SaveClonePNGHR.cpp
        ATCG/SaveClonePNGS.cpp
        ATCG/SavePNGHR.cpp
        ATCG/SavePNGS.cpp
        ATCG/SavePNGSingleCell.cpp
        ATCG/SaveVisualArrayL1.cpp
        ATCG/SaveVisualArrayL2.cpp
        ATCG/SaveVisualArrayL3.cpp
        ATCG/cell_store.cpp
        ATCG/cell_trace.cpp
        ATCG/cell_motion.cpp
        ATCG/color_space.cpp
        ATCG/cell_type_transform.cpp
        ATCG/common/density_growth_rule.cpp
        ATCG/death_judgement.cpp
        ATCG/deltah_calculation.cpp
        ATCG/deltah_recalculation.cpp
        ATCG/density_calculation.cpp
        ATCG/density_growth_rate_calculation_1.cpp
        ATCG/density_dependent_growth.cpp
        ATCG/division.cpp
        ATCG/free_living_division.cpp
        ATCG/free_living_division_single_thread.cpp
        ATCG/free_living_growth.cpp
        ATCG/free_living_growth_recovery.cpp
        ATCG/inner_count.cpp
        ATCG/inner_count_low_density.cpp
        ATCG/inner_initiation_array.cpp
        ATCG/int_grid.cpp
        ATCG/int_matrix.cpp
        ATCG/main.cpp
        ATCG/Low_density_initial_growth.cpp
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
        ATCG/save_data.cpp
        ATCG/save_data_free_living.cpp
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
endif()

if(ATCG_BUILD_3D)
    add_library(atcg3d_core
        ATCG/common/density_growth_rule.cpp
        ATCG3D/config/model_config.cpp
        ATCG3D/core/cell_store.cpp
        ATCG3D/engine/simulation.cpp
        ATCG3D/io/output_manager.cpp
        ATCG3D/io/preview_sampler.cpp
        ATCG3D/io/run_manifest.cpp
        ATCG3D/rules/density.cpp
        ATCG3D/rules/initialization.cpp
        ATCG3D/rules/lifecycle.cpp
        ATCG3D/rules/migration.cpp
        ATCG3D/space/chunk_grid.cpp
        ATCG3D/space/density_index.cpp
    )
    target_include_directories(atcg3d_core PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/ATCG3D
        ${CMAKE_CURRENT_SOURCE_DIR}/ATCG
    )
    target_link_libraries(atcg3d_core PUBLIC atcg_openmp)
    target_compile_definitions(atcg3d_core PUBLIC ATCG3D_SCHEMA_VERSION=1)

    add_executable(atcg3d ATCG3D/app/main.cpp)
    target_link_libraries(atcg3d PRIVATE atcg3d_core)

    add_executable(atcg3d_scale_benchmark benchmarks/atcg3d_scale.cpp)
    target_link_libraries(atcg3d_scale_benchmark PRIVATE atcg3d_core)
    add_executable(atcg3d_event_benchmark benchmarks/atcg3d_events.cpp)
    target_link_libraries(atcg3d_event_benchmark PRIVATE atcg3d_core)

    file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/configs)
    configure_file(
        ${CMAKE_CURRENT_SOURCE_DIR}/configs/atcg3d_legacy_like_v1.cfg
        ${CMAKE_CURRENT_BINARY_DIR}/configs/atcg3d_legacy_like_v1.cfg
        COPYONLY
    )

    if(ATCG3D_ENABLE_VTKHDF)
        find_package(VTK REQUIRED COMPONENTS CommonCore CommonDataModel IOHDF)
        target_sources(atcg3d_core PRIVATE ATCG3D/io/vtkhdf_writer.cpp)
        target_link_libraries(atcg3d_core PUBLIC
            VTK::CommonCore
            VTK::CommonDataModel
            VTK::IOHDF
        )
        target_compile_definitions(atcg3d_core PUBLIC ATCG3D_HAS_VTKHDF=1)
        vtk_module_autoinit(
            TARGETS atcg3d_core atcg3d
            MODULES VTK::CommonCore VTK::CommonDataModel VTK::IOHDF
        )
        message(STATUS "VTK-HDF visualization output enabled")
    else()
        target_sources(atcg3d_core PRIVATE ATCG3D/io/vtkhdf_writer_unavailable.cpp)
    endif()

    if(ATCG3D_ENABLE_HDF5_CHECKPOINT)
        find_package(HDF5 REQUIRED COMPONENTS CXX)
        target_sources(atcg3d_core PRIVATE ATCG3D/io/checkpoint_hdf5.cpp)
        target_include_directories(atcg3d_core PUBLIC ${HDF5_INCLUDE_DIRS})
        target_link_libraries(atcg3d_core PUBLIC ${HDF5_CXX_LIBRARIES})
        target_compile_definitions(atcg3d_core PUBLIC ATCG3D_HAS_HDF5_CHECKPOINT=1)
        message(STATUS "HDF5 checkpoint support enabled")
    endif()
endif()

include(CTest)

if(BUILD_TESTING)
    find_package(Python3 QUIET COMPONENTS Interpreter)
    add_executable(stateless_rng_test
        ATCG/stateless_rng.cpp
        tests/stateless_rng_test.cpp
    )
    target_include_directories(stateless_rng_test PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/ATCG)
    add_test(NAME stateless_rng_test COMMAND stateless_rng_test)

    add_executable(density_growth_rate_test
        ATCG/cell_store.cpp
        ATCG/common/density_growth_rule.cpp
        ATCG/density_growth_rate_calculation_1.cpp
        ATCG/stateless_rng.cpp
        ATCG/visual_range.cpp
        tests/density_growth_rate_test.cpp
    )
    target_include_directories(density_growth_rate_test PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/ATCG)
    add_test(NAME density_growth_rate_test COMMAND density_growth_rate_test)

    if(ATCG_BUILD_3D)
        foreach(test_name IN ITEMS geometry config storage density migration lifecycle simulation io_core)
            add_executable(atcg3d_${test_name}_test tests/3d/${test_name}_test.cpp)
            target_link_libraries(atcg3d_${test_name}_test PRIVATE atcg3d_core)
            add_test(NAME atcg3d_${test_name}_test COMMAND atcg3d_${test_name}_test)
        endforeach()
        if(ATCG3D_ENABLE_HDF5_CHECKPOINT)
            add_executable(atcg3d_checkpoint_test tests/3d/checkpoint_test.cpp)
            target_link_libraries(atcg3d_checkpoint_test PRIVATE atcg3d_core)
            add_test(NAME atcg3d_checkpoint_test COMMAND atcg3d_checkpoint_test)
        endif()
        if(ATCG3D_ENABLE_VTKHDF)
            add_executable(atcg3d_vtkhdf_test tests/3d/vtkhdf_test.cpp)
            target_link_libraries(atcg3d_vtkhdf_test PRIVATE atcg3d_core)
            vtk_module_autoinit(
                TARGETS atcg3d_vtkhdf_test
                MODULES VTK::CommonCore VTK::CommonDataModel VTK::IOHDF
            )
            add_test(NAME atcg3d_vtkhdf_test COMMAND atcg3d_vtkhdf_test)
        endif()
        if(Python3_Interpreter_FOUND)
            add_test(
                NAME atcg3d_viewer_controller_test
                COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/3d/viewer_controller_test.py
            )
            set_tests_properties(atcg3d_viewer_controller_test PROPERTIES
                WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
            )
        endif()
        if(ATCG3D_ENABLE_LARGE_TESTS)
            add_test(
                NAME atcg3d_ten_million_core_benchmark
                COMMAND atcg3d_scale_benchmark
                    --cells 10000000
                    --skip-vtk
                    --skip-checkpoint
                    --directory ${CMAKE_CURRENT_BINARY_DIR}/atcg3d_10m_core
                    --result ${CMAKE_CURRENT_BINARY_DIR}/atcg3d_10m_core/result.json
            )
            set_tests_properties(atcg3d_ten_million_core_benchmark PROPERTIES LABELS large)
        endif()
    endif()
endif()
