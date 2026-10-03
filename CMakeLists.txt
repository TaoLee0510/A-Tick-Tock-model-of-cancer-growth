cmake_minimum_required(VERSION 3.20)

project(ATCG LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS ON)

option(ATCG_ENABLE_LONG_VALIDATION "Register opt-in long ensembles" OFF)
option(ATCG_BUILD_LEGACY_2D "Build the existing 2D ATCG executable" ON)
option(ATCG_BUILD_3D "Build the sparse event-driven 3D executable" ON)
option(ATCG_BUILD_3D_NUTRIENT "Build the nutrient-coupled 3D executable" ON)
option(ATCG_BUILD_3D_CONTINUUM "Build the population-density continuum executable" ON)
option(ATCG_BUILD_3D_STRUCTURED_PDE "Build the ABM-aligned structured PDE executable" ON)
option(ATCG3D_ENABLE_VTKHDF "Build VTK-HDF visualization output" OFF)
option(ATCG3D_ENABLE_HDF5_CHECKPOINT "Build HDF5 checkpoint support" OFF)
option(ATCG3D_ENABLE_LARGE_TESTS "Enable opt-in million/ten-million-cell tests" OFF)

if(APPLE)
    find_program(ATCG_HOMEBREW_EXECUTABLE brew)
    if(ATCG_HOMEBREW_EXECUTABLE)
        execute_process(COMMAND ${ATCG_HOMEBREW_EXECUTABLE} --prefix
            OUTPUT_VARIABLE ATCG_HOMEBREW_PREFIX
            OUTPUT_STRIP_TRAILING_WHITESPACE
            RESULT_VARIABLE ATCG_HOMEBREW_RESULT
            ERROR_QUIET)
        if(ATCG_HOMEBREW_RESULT EQUAL 0 AND ATCG_HOMEBREW_PREFIX)
            list(APPEND CMAKE_PREFIX_PATH
                "${ATCG_HOMEBREW_PREFIX}" "${ATCG_HOMEBREW_PREFIX}/opt/libomp")
        endif()
    endif()
endif()

if(ATCG_BUILD_LEGACY_2D OR ATCG_BUILD_3D)
    if(APPLE AND CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        find_path(OPENMP_INCLUDE_DIR
            NAMES omp.h
            REQUIRED
        )
        find_library(OPENMP_LIBRARY
            NAMES omp
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
        REQUIRED
    )

    find_path(PNGWRITER_INCLUDE_DIR
        NAMES pngwriter.h
        REQUIRED
    )

    find_library(PNGWRITER_LIBRARY
        NAMES pngwriter
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
    find_package(yaml-cpp 0.7 REQUIRED CONFIG)
    if(TARGET yaml-cpp::yaml-cpp)
        set(ATCG3D_YAML_TARGET yaml-cpp::yaml-cpp)
    elseif(TARGET yaml-cpp)
        set(ATCG3D_YAML_TARGET yaml-cpp)
    else()
        message(FATAL_ERROR "yaml-cpp did not export a supported CMake target")
    endif()

    add_library(atcg3d_core
        ATCG/common/density_growth_rule.cpp
        ATCG3D/config/model_config.cpp
        ATCG3D/config/yaml_loader.cpp
        ATCG3D/app/run_controller.cpp
        ATCG3D/core/cell_store.cpp
        ATCG3D/engine/parallelism.cpp
        ATCG3D/engine/simulation.cpp
        ATCG3D/io/output_manager.cpp
        ATCG3D/io/preview_sampler.cpp
        ATCG3D/io/run_manifest.cpp
        ATCG3D/rules/density.cpp
        ATCG3D/rules/initialization.cpp
        ATCG3D/rules/initial_rates.cpp
        ATCG3D/rules/lifecycle.cpp
        ATCG3D/rules/migration.cpp
        ATCG3D/space/chunk_grid.cpp
        ATCG3D/space/density_index.cpp
        ATCG3D/space/quantized_box_count_index.cpp
        ATCG3D/vasculature/angiogenesis_process.cpp
        ATCG3D/vasculature/geometry.cpp
        ATCG3D/vasculature/growth_rules.cpp
        ATCG3D/vasculature/influence_field.cpp
        ATCG3D/vasculature/lesion_index.cpp
        ATCG3D/vasculature/surface_index.cpp
        ATCG3D/vasculature/vessel_grid.cpp
        ATCG3D/vasculature/vessel_store.cpp
    )
    target_include_directories(atcg3d_core PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/ATCG3D
        ${CMAKE_CURRENT_SOURCE_DIR}/ATCG
    )
    target_link_libraries(atcg3d_core PUBLIC atcg_openmp PRIVATE ${ATCG3D_YAML_TARGET})
    target_compile_definitions(atcg3d_core PUBLIC ATCG3D_SCHEMA_VERSION=3)

    add_executable(atcg3d ATCG3D/app/main.cpp)
    target_link_libraries(atcg3d PRIVATE atcg3d_core)

    if(ATCG_BUILD_3D_NUTRIENT)
        add_library(atcg3d_nutrient_core
            ATCG3D_Nutrient/config/nutrient_config.cpp
            ATCG3D_Nutrient/field/nutrient_field.cpp
            ATCG3D_Nutrient/io/nutrient_output.cpp
        )
        target_include_directories(atcg3d_nutrient_core PUBLIC
            ${CMAKE_CURRENT_SOURCE_DIR}/ATCG3D_Nutrient
        )
        target_link_libraries(atcg3d_nutrient_core
            PUBLIC atcg3d_core
            PRIVATE ${ATCG3D_YAML_TARGET})

        add_executable(atcg3d_nutrient ATCG3D_Nutrient/app/main.cpp)
        target_link_libraries(atcg3d_nutrient PRIVATE atcg3d_nutrient_core)

        foreach(nutrient_config_name IN ITEMS
                nutrient_smoke_v1.yaml
                nutrient_legacy_v1.yaml
                nutrient_smoke_per_cell_r12_guided_v2.yaml
                nutrient_legacy_per_cell_r12_guided_v2.yaml)
            configure_file(
                ${CMAKE_CURRENT_SOURCE_DIR}/ATCG3D_Nutrient/config/${nutrient_config_name}
                ${CMAKE_CURRENT_BINARY_DIR}/configs/${nutrient_config_name}
                COPYONLY
            )
        endforeach()
    endif()

    if(ATCG_BUILD_3D_CONTINUUM)
        add_library(atcg3d_continuum_core
            ATCG3D_Continuum/model/angiogenesis_field.cpp
            ATCG3D_Continuum/model/shared_angiogenesis.cpp
            ATCG3D_Continuum/model/sector_mean_field.cpp
            ATCG3D_Continuum/config/continuum_config.cpp
            ATCG3D_Continuum/model/moving_tumor_front.cpp
            ATCG3D_Continuum/model/continuum_model.cpp
            ATCG3D_Continuum/io/continuum_output.cpp
            ATCG3D_Continuum/io/pde_vtkhdf_writer.cpp
        )
        target_include_directories(atcg3d_continuum_core PUBLIC
            ${CMAKE_CURRENT_SOURCE_DIR}/ATCG3D_Continuum
        )
        target_link_libraries(atcg3d_continuum_core
            PUBLIC atcg3d_core
            PRIVATE ${ATCG3D_YAML_TARGET})

        add_executable(atcg3d_continuum ATCG3D_Continuum/app/main.cpp)
        target_link_libraries(atcg3d_continuum PRIVATE atcg3d_continuum_core)
        add_executable(atcg3d_angiogenesis_benchmark ATCG3D_Continuum/app/angiogenesis_benchmark.cpp)
        target_link_libraries(atcg3d_angiogenesis_benchmark PRIVATE atcg3d_continuum_core)

        foreach(continuum_config_name IN ITEMS
                continuum_smoke_v1.yaml
                continuum_legacy_v1.yaml)
            configure_file(
                ${CMAKE_CURRENT_SOURCE_DIR}/ATCG3D_Continuum/config/${continuum_config_name}
                ${CMAKE_CURRENT_BINARY_DIR}/configs/${continuum_config_name}
                COPYONLY
            )
        endforeach()
    endif()

    if(ATCG_BUILD_3D_STRUCTURED_PDE)
        if(NOT ATCG_BUILD_3D_CONTINUUM)
            message(FATAL_ERROR
                "ATCG_BUILD_3D_STRUCTURED_PDE requires ATCG_BUILD_3D_CONTINUUM")
        endif()
        add_library(atcg3d_structured_pde_core
            ATCG3D_StructuredPDE/config/structured_config.cpp
            ATCG3D_StructuredPDE/model/structured_pde_model.cpp
            ATCG3D_StructuredPDE/model/feasible_jump.cpp
            ATCG3D_StructuredPDE/model/division_renewal.cpp
            ATCG3D_StructuredPDE/model/beta_duration.cpp
            ATCG3D_StructuredPDE/io/structured_pde_output.cpp
            ATCG3D_StructuredPDE/io/structured_vtkhdf_writer.cpp
        )
        target_include_directories(atcg3d_structured_pde_core PUBLIC
            ${CMAKE_CURRENT_SOURCE_DIR}/ATCG3D_StructuredPDE
        )
        target_link_libraries(atcg3d_structured_pde_core
            PUBLIC atcg3d_continuum_core
            PRIVATE ${ATCG3D_YAML_TARGET})

        add_executable(atcg3d_structured_pde
            ATCG3D_StructuredPDE/app/main.cpp)
        add_library(atcg3d_shared_rules_core
            ATCG3D_SharedRules/model/shared_resource_environment.cpp)
        target_include_directories(atcg3d_shared_rules_core PUBLIC
            ${CMAKE_CURRENT_SOURCE_DIR}/ATCG3D_SharedRules)
        target_link_libraries(atcg3d_shared_rules_core PUBLIC atcg3d_structured_pde_core)
        target_link_libraries(atcg3d_structured_pde PRIVATE atcg3d_shared_rules_core)
        add_executable(atcg3d_shared_abm ATCG3D_SharedRules/app/main.cpp)
        target_link_libraries(atcg3d_shared_abm PRIVATE atcg3d_shared_rules_core)
        add_executable(atcg3d_shared_3d_benchmark benchmarks/atcg3d_shared_3d.cpp)
        target_link_libraries(atcg3d_shared_3d_benchmark PRIVATE atcg3d_shared_rules_core)

        add_library(atcg3d_hybrid_core
            ATCG3D_Hybrid/model/hybrid_model.cpp
            ATCG3D_Hybrid/model/hybrid_invasion.cpp
            ATCG3D_Hybrid/model/hybrid_conversion.cpp)
        target_include_directories(atcg3d_hybrid_core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/ATCG3D_Hybrid)
        target_link_libraries(atcg3d_hybrid_core PUBLIC atcg3d_shared_rules_core PRIVATE ${ATCG3D_YAML_TARGET})
        add_executable(atcg3d_hybrid ATCG3D_Hybrid/app/main.cpp)
        target_link_libraries(atcg3d_hybrid PRIVATE atcg3d_hybrid_core)

        add_library(atcg3d_ode_core ATCG3D_ODE/ode/ode_model.cpp)
        target_include_directories(atcg3d_ode_core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/ATCG3D_ODE)
        target_link_libraries(atcg3d_ode_core PUBLIC atcg3d_structured_pde_core PRIVATE ${ATCG3D_YAML_TARGET})
        add_executable(atcg3d_ode ATCG3D_ODE/app/main.cpp)
        target_link_libraries(atcg3d_ode PRIVATE atcg3d_ode_core)

        add_executable(atcg_sim apps/atcg_sim.cpp)
        target_link_libraries(atcg_sim PRIVATE ${ATCG3D_YAML_TARGET})
        add_executable(atcg_config_migrate apps/config_migrate.cpp ATCG3D_Nutrient/config/nutrient_config.cpp)
        target_include_directories(atcg_config_migrate PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/ATCG3D_Nutrient)
        target_link_libraries(atcg_config_migrate PRIVATE atcg3d_hybrid_core atcg3d_ode_core ${ATCG3D_YAML_TARGET})
        install(TARGETS atcg_sim atcg_config_migrate atcg3d atcg3d_continuum atcg3d_structured_pde atcg3d_shared_abm atcg3d_ode atcg3d_hybrid RUNTIME DESTINATION bin)
        if(TARGET atcg3d_nutrient)
            install(TARGETS atcg3d_nutrient RUNTIME DESTINATION bin)
        endif()
        install(DIRECTORY docs DESTINATION share/atcg)
        install(DIRECTORY configs ATCG3D_Nutrient ATCG3D_Continuum ATCG3D_StructuredPDE ATCG3D_StructuredPDE_NutrientChemotaxis ATCG3D_SharedRules ATCG3D_ODE ATCG3D_Hybrid DESTINATION share/atcg FILES_MATCHING PATTERN "*.yaml")

        foreach(structured_config_name IN ITEMS
                continuum_legacy_2d_2000_r20_aligned_v1.yaml
                structured_legacy_2d_2000_r20_v1.yaml
                continuum_smoke_r20_v1.yaml
                structured_smoke_r20_v1.yaml
                continuum_benchmark_2d_2000_24h_r20_v1.yaml
                structured_benchmark_2d_2000_24h_r20_v1.yaml
                continuum_legacy_2d_2000_r20_resource_guided_v2.yaml
                structured_legacy_2d_2000_r20_resource_guided_v2.yaml
                continuum_smoke_r20_resource_guided_v2.yaml
                structured_smoke_r20_resource_guided_v2.yaml
                continuum_benchmark_2d_2000_24h_resource_guided_v2.yaml
                structured_benchmark_2d_2000_24h_resource_guided_v2.yaml
                continuum_legacy_2d_2000_r20_exchange_v3.yaml
                structured_legacy_2d_2000_r20_exchange_v3.yaml
                continuum_smoke_r20_exchange_v3.yaml
                structured_smoke_r20_exchange_v3.yaml)
            configure_file(
                ${CMAKE_CURRENT_SOURCE_DIR}/ATCG3D_StructuredPDE/config/${structured_config_name}
                ${CMAKE_CURRENT_BINARY_DIR}/configs/${structured_config_name}
                COPYONLY
            )
        endforeach()
    endif()

    add_executable(atcg3d_scale_benchmark benchmarks/atcg3d_scale.cpp)
    target_link_libraries(atcg3d_scale_benchmark PRIVATE atcg3d_core)
    add_executable(atcg3d_event_benchmark benchmarks/atcg3d_events.cpp)
    target_link_libraries(atcg3d_event_benchmark PRIVATE atcg3d_core)
    add_executable(atcg3d_initialization_benchmark benchmarks/atcg3d_initialization.cpp)
    target_link_libraries(atcg3d_initialization_benchmark PRIVATE atcg3d_core)
    add_executable(atcg3d_angiogenesis_scale_benchmark
        benchmarks/atcg3d_angiogenesis_scale.cpp)
    target_link_libraries(atcg3d_angiogenesis_scale_benchmark PRIVATE atcg3d_core)

    file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/configs)
    # Remove stale schema-v2-named copies from an existing build tree so a
    # renamed profile cannot be selected accidentally after reconfiguration.
    file(REMOVE
        ${CMAKE_CURRENT_BINARY_DIR}/configs/atcg3d_legacy_2d_mapped_v2.yaml
        ${CMAKE_CURRENT_BINARY_DIR}/configs/atcg3d_smoke_test_v2.yaml
    )
    foreach(config_name IN ITEMS
            atcg3d_legacy_2d_mapped_v3.yaml
            atcg3d_smoke_test_v3.yaml
            single_r_stage0_2160h_density_vascular_v5.yaml
            single_r_stage0_2160h_seed1.yaml
            atcg2d_legacy_native_r20_aligned_v3.yaml)
        configure_file(
            ${CMAKE_CURRENT_SOURCE_DIR}/configs/${config_name}
            ${CMAKE_CURRENT_BINARY_DIR}/configs/${config_name}
            COPYONLY
        )
    endforeach()

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
            TARGETS
                atcg3d_core
                atcg3d
                atcg3d_scale_benchmark
                atcg3d_event_benchmark
                atcg3d_initialization_benchmark
                atcg3d_angiogenesis_scale_benchmark
            MODULES VTK::CommonCore VTK::CommonDataModel VTK::IOHDF
        )
        if(ATCG_BUILD_3D_NUTRIENT)
            vtk_module_autoinit(
                TARGETS atcg3d_nutrient
                MODULES VTK::CommonCore VTK::CommonDataModel VTK::IOHDF
            )
        endif()
        if(ATCG_BUILD_3D_CONTINUUM)
            vtk_module_autoinit(
                TARGETS atcg3d_continuum
                MODULES VTK::CommonCore VTK::CommonDataModel VTK::IOHDF
            )
        endif()
        if(ATCG_BUILD_3D_STRUCTURED_PDE)
            vtk_module_autoinit(
                TARGETS atcg3d_structured_pde
                MODULES VTK::CommonCore VTK::CommonDataModel VTK::IOHDF
            )
        endif()
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
    if(ATCG_BUILD_3D_CONTINUUM)
        add_executable(atcg3d_sector_mean_test tests/3d/sector_mean_test.cpp)
        target_link_libraries(atcg3d_sector_mean_test PRIVATE atcg3d_continuum_core)
        target_compile_options(atcg3d_sector_mean_test PRIVATE -UNDEBUG)
        add_test(NAME atcg3d_sector_mean_test COMMAND atcg3d_sector_mean_test)
        set_tests_properties(atcg3d_sector_mean_test PROPERTIES TIMEOUT 600)
    endif()
    if(ATCG_BUILD_3D_STRUCTURED_PDE)
        add_executable(atcg3d_three_dimensional_shared_test tests/3d/three_dimensional_shared_test.cpp)
        target_link_libraries(atcg3d_three_dimensional_shared_test PRIVATE atcg3d_shared_rules_core)
        target_compile_options(atcg3d_three_dimensional_shared_test PRIVATE -UNDEBUG)
        target_compile_definitions(atcg3d_three_dimensional_shared_test PRIVATE ATCG_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")
        add_test(NAME atcg3d_three_dimensional_shared_test COMMAND atcg3d_three_dimensional_shared_test)
        set_tests_properties(atcg3d_three_dimensional_shared_test PROPERTIES TIMEOUT 600)
    endif()
    find_package(Python3 QUIET COMPONENTS Interpreter)
    if(Python3_Interpreter_FOUND AND ATCG_BUILD_3D_STRUCTURED_PDE AND ATCG3D_ENABLE_LARGE_TESTS)
        add_test(NAME atcg3d_shared_3d_benchmark
            COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/scripts/benchmark_shared_3d.py
                --exe $<TARGET_FILE:atcg3d_shared_3d_benchmark>
                --output ${CMAKE_CURRENT_BINARY_DIR}/benchmark-shared-3d.json --repeat 1)
        set_tests_properties(atcg3d_shared_3d_benchmark PROPERTIES LABELS "benchmark;benchmark_3d" TIMEOUT 3600)
    endif()
    add_executable(stateless_rng_test
        ATCG/stateless_rng.cpp
        tests/stateless_rng_test.cpp
    )
    target_include_directories(stateless_rng_test PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/ATCG)
    target_compile_options(stateless_rng_test PRIVATE -UNDEBUG)
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
    target_compile_options(density_growth_rate_test PRIVATE -UNDEBUG)
    add_test(NAME density_growth_rate_test COMMAND density_growth_rate_test)

    if(ATCG_BUILD_3D)
        foreach(test_name IN ITEMS
                geometry config storage density migration lifecycle initialization initial_rates simulation io_core
                angiogenesis_geometry angiogenesis_process angiogenesis_direction parallelism
                angiogenesis_simulation lesion_index vascular_density run_controller)
            add_executable(atcg3d_${test_name}_test tests/3d/${test_name}_test.cpp)
            target_link_libraries(atcg3d_${test_name}_test PRIVATE atcg3d_core)
            target_compile_options(atcg3d_${test_name}_test PRIVATE -UNDEBUG)
            add_test(NAME atcg3d_${test_name}_test COMMAND atcg3d_${test_name}_test)
        endforeach()
        if(ATCG_BUILD_3D_NUTRIENT)
            add_executable(atcg3d_nutrient_test
                tests/3d/nutrient_test.cpp)
            target_link_libraries(atcg3d_nutrient_test PRIVATE
                atcg3d_nutrient_core)
            target_compile_options(atcg3d_nutrient_test PRIVATE -UNDEBUG)
            target_compile_definitions(atcg3d_nutrient_test PRIVATE
                ATCG_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")
            add_test(NAME atcg3d_nutrient_test COMMAND atcg3d_nutrient_test)
        endif()
        if(ATCG_BUILD_3D_CONTINUUM)
            add_executable(atcg3d_continuum_test
                tests/3d/continuum_test.cpp)
            target_link_libraries(atcg3d_continuum_test PRIVATE
                atcg3d_continuum_core)
            target_compile_options(atcg3d_continuum_test PRIVATE -UNDEBUG)
            target_compile_definitions(atcg3d_continuum_test PRIVATE
                ATCG_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")
            add_test(NAME atcg3d_continuum_test COMMAND atcg3d_continuum_test)
        endif()
        if(ATCG_BUILD_3D_STRUCTURED_PDE)
            add_executable(atcg3d_shared_rules_test tests/3d/shared_rules_test.cpp)
            target_link_libraries(atcg3d_shared_rules_test PRIVATE atcg3d_shared_rules_core)
            target_compile_options(atcg3d_shared_rules_test PRIVATE -UNDEBUG)
            target_compile_definitions(atcg3d_shared_rules_test PRIVATE
                ATCG_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")
            add_test(NAME atcg3d_shared_rules_test COMMAND atcg3d_shared_rules_test)
            add_executable(atcg3d_division_renewal_test tests/3d/division_renewal_test.cpp)
            target_link_libraries(atcg3d_division_renewal_test PRIVATE atcg3d_structured_pde_core)
            target_compile_options(atcg3d_division_renewal_test PRIVATE -UNDEBUG)
            target_compile_definitions(atcg3d_division_renewal_test PRIVATE ATCG_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")
            add_test(NAME atcg3d_division_renewal_test COMMAND atcg3d_division_renewal_test)
            add_executable(atcg3d_activation_distribution_test tests/3d/activation_distribution_test.cpp)
            target_link_libraries(atcg3d_activation_distribution_test PRIVATE atcg3d_structured_pde_core)
            target_compile_options(atcg3d_activation_distribution_test PRIVATE -UNDEBUG)
            target_compile_definitions(atcg3d_activation_distribution_test PRIVATE ATCG_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")
            add_test(NAME atcg3d_activation_distribution_test COMMAND atcg3d_activation_distribution_test)
            add_executable(atcg3d_product_cli_test tests/3d/product_cli_test.cpp)
            target_link_libraries(atcg3d_product_cli_test PRIVATE atcg3d_hybrid_core atcg3d_ode_core ${ATCG3D_YAML_TARGET})
            target_compile_options(atcg3d_product_cli_test PRIVATE -UNDEBUG)
            target_compile_definitions(atcg3d_product_cli_test PRIVATE ATCG_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}" ATCG_BINARY_DIR="${CMAKE_CURRENT_BINARY_DIR}")
            add_dependencies(atcg3d_product_cli_test atcg_sim atcg_config_migrate atcg3d_structured_pde atcg3d_shared_abm atcg3d_ode atcg3d_hybrid)
            add_test(NAME atcg3d_product_cli_test COMMAND atcg3d_product_cli_test)
            if(ATCG3D_ENABLE_HDF5_CHECKPOINT)
                add_executable(atcg3d_pde_vtkhdf_test tests/3d/pde_vtkhdf_test.cpp)
                target_link_libraries(atcg3d_pde_vtkhdf_test PRIVATE atcg3d_structured_pde_core)
                target_compile_options(atcg3d_pde_vtkhdf_test PRIVATE -UNDEBUG)
                target_compile_definitions(atcg3d_pde_vtkhdf_test PRIVATE ATCG_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")
                add_test(NAME atcg3d_pde_vtkhdf_test COMMAND atcg3d_pde_vtkhdf_test)
            endif()
            add_executable(atcg3d_sparse_pde_test tests/3d/sparse_pde_test.cpp)
            target_link_libraries(atcg3d_sparse_pde_test PRIVATE atcg3d_structured_pde_core)
            target_compile_options(atcg3d_sparse_pde_test PRIVATE -UNDEBUG)
            target_compile_definitions(atcg3d_sparse_pde_test PRIVATE ATCG_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")
            add_test(NAME atcg3d_sparse_pde_test COMMAND atcg3d_sparse_pde_test)
            if(ATCG3D_ENABLE_LARGE_TESTS)
                add_test(NAME atcg3d_sparse_pde_10000_memory COMMAND atcg3d_sparse_pde_test --large)
                set_tests_properties(atcg3d_sparse_pde_10000_memory PROPERTIES LABELS benchmark TIMEOUT 120)
            endif()
            add_executable(atcg3d_hybrid_test tests/3d/hybrid_test.cpp)
            target_link_libraries(atcg3d_hybrid_test PRIVATE atcg3d_hybrid_core)
            target_compile_options(atcg3d_hybrid_test PRIVATE -UNDEBUG)
            target_compile_definitions(atcg3d_hybrid_test PRIVATE ATCG_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")
            add_test(NAME atcg3d_hybrid_test COMMAND atcg3d_hybrid_test)
            add_executable(atcg3d_hybrid_distribution_test tests/3d/hybrid_distribution_test.cpp)
            target_link_libraries(atcg3d_hybrid_distribution_test PRIVATE atcg3d_hybrid_core)
            target_compile_options(atcg3d_hybrid_distribution_test PRIVATE -UNDEBUG)
            target_compile_definitions(atcg3d_hybrid_distribution_test PRIVATE ATCG_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")
            add_test(NAME atcg3d_hybrid_distribution_test COMMAND atcg3d_hybrid_distribution_test)
            add_executable(atcg3d_hybrid_volume_test tests/3d/hybrid_volume_test.cpp)
            target_link_libraries(atcg3d_hybrid_volume_test PRIVATE atcg3d_hybrid_core)
            target_compile_options(atcg3d_hybrid_volume_test PRIVATE -UNDEBUG)
            target_compile_definitions(atcg3d_hybrid_volume_test PRIVATE ATCG_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")
            add_test(NAME atcg3d_hybrid_volume_test COMMAND atcg3d_hybrid_volume_test)
            add_executable(atcg3d_hybrid_invasion_test tests/3d/hybrid_invasion_test.cpp)
            target_link_libraries(atcg3d_hybrid_invasion_test PRIVATE atcg3d_hybrid_core)
            target_compile_options(atcg3d_hybrid_invasion_test PRIVATE -UNDEBUG)
            target_compile_definitions(atcg3d_hybrid_invasion_test PRIVATE ATCG_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")
            add_test(NAME atcg3d_hybrid_invasion_test COMMAND atcg3d_hybrid_invasion_test)
            add_executable(atcg3d_vascular_removal_test tests/3d/vascular_removal_test.cpp)
            target_link_libraries(atcg3d_vascular_removal_test PRIVATE atcg3d_structured_pde_core)
            target_compile_options(atcg3d_vascular_removal_test PRIVATE -UNDEBUG)
            target_compile_definitions(atcg3d_vascular_removal_test PRIVATE ATCG_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")
            add_test(NAME atcg3d_vascular_removal_test COMMAND atcg3d_vascular_removal_test)
            add_executable(atcg3d_ode_test tests/3d/ode_test.cpp)
            target_link_libraries(atcg3d_ode_test PRIVATE atcg3d_ode_core)
            target_compile_options(atcg3d_ode_test PRIVATE -UNDEBUG)
            target_compile_definitions(atcg3d_ode_test PRIVATE ATCG_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")
            add_test(NAME atcg3d_ode_test COMMAND atcg3d_ode_test)
            add_executable(atcg3d_vascular_field_test tests/3d/vascular_field_test.cpp)
            target_link_libraries(atcg3d_vascular_field_test PRIVATE atcg3d_shared_rules_core)
            target_compile_options(atcg3d_vascular_field_test PRIVATE -UNDEBUG)
            target_compile_definitions(atcg3d_vascular_field_test PRIVATE ATCG_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")
            add_test(NAME atcg3d_vascular_field_test COMMAND atcg3d_vascular_field_test)
            add_executable(atcg3d_shared_angiogenesis_test tests/3d/shared_angiogenesis_test.cpp)
            target_link_libraries(atcg3d_shared_angiogenesis_test PRIVATE atcg3d_shared_rules_core)
            target_compile_options(atcg3d_shared_angiogenesis_test PRIVATE -UNDEBUG)
            add_test(NAME atcg3d_shared_angiogenesis_test COMMAND atcg3d_shared_angiogenesis_test)
            add_executable(atcg3d_shared_vascular_coupling_test tests/3d/shared_vascular_coupling_test.cpp)
            target_link_libraries(atcg3d_shared_vascular_coupling_test PRIVATE atcg3d_shared_rules_core)
            target_compile_options(atcg3d_shared_vascular_coupling_test PRIVATE -UNDEBUG)
            target_compile_definitions(atcg3d_shared_vascular_coupling_test PRIVATE ATCG_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")
            add_test(NAME atcg3d_shared_vascular_coupling_test COMMAND atcg3d_shared_vascular_coupling_test)
            add_executable(atcg3d_operator_contract_test tests/3d/operator_contract_test.cpp)
            target_link_libraries(atcg3d_operator_contract_test PRIVATE atcg3d_shared_rules_core)
            target_compile_options(atcg3d_operator_contract_test PRIVATE -UNDEBUG)
            target_compile_definitions(atcg3d_operator_contract_test PRIVATE ATCG_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")
            add_test(NAME atcg3d_operator_contract_test COMMAND atcg3d_operator_contract_test)
            if(Python3_Interpreter_FOUND)
                add_test(NAME atcg3d_equivalence_stats_test
                    COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/3d/equivalence_stats_test.py)
                add_test(NAME atcg3d_hybrid_invasion_stats_test
                    COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/3d/hybrid_invasion_stats_test.py)
                add_test(NAME atcg3d_hybrid_invasion_trace_test
                    COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/3d/hybrid_invasion_trace_test.py
                        --exe $<TARGET_FILE:atcg3d_hybrid>
                        --native-exe $<TARGET_FILE:atcg3d_shared_abm>
                        --root ${CMAKE_CURRENT_SOURCE_DIR})
                add_test(NAME atcg3d_shared_angiogenesis_stats_test
                    COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/3d/shared_angiogenesis_stats_test.py)
                add_test(NAME atcg3d_validation_trace_test
                    COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/3d/validation_trace_test.py
                        --exe $<TARGET_FILE:atcg3d_shared_abm> --root ${CMAKE_CURRENT_SOURCE_DIR})
                if(ATCG3D_ENABLE_HDF5_CHECKPOINT)
                    add_test(NAME atcg3d_validation_trace_checkpoint_test
                        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/3d/validation_trace_test.py
                            --exe $<TARGET_FILE:atcg3d_shared_abm> --root ${CMAKE_CURRENT_SOURCE_DIR} --checkpoint)
                endif()
                add_test(NAME atcg3d_statistical_regular_validation
                    COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/scripts/validate_abm_vs_pde.py
                        --exe $<TARGET_FILE:atcg3d_shared_abm>
                        --config ${CMAKE_CURRENT_SOURCE_DIR}/ATCG3D_SharedRules/config/regular_cycle_spatial_birth_v14.yaml
                        --seeds 16 --require-interior
                        --output ${CMAKE_CURRENT_BINARY_DIR}/validation-statistical-spatial-birth)
                set_tests_properties(atcg3d_statistical_regular_validation PROPERTIES LABELS validation TIMEOUT 600)
                add_test(NAME atcg3d_statistical_active_r20_validation
                    COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/scripts/validate_abm_vs_pde.py
                        --exe $<TARGET_FILE:atcg3d_shared_abm>
                        --config ${CMAKE_CURRENT_SOURCE_DIR}/ATCG3D_SharedRules/config/active_r20_v14.yaml
                        --seeds 16 --require-interior --minimum-active-fraction 0.1
                        --output ${CMAKE_CURRENT_BINARY_DIR}/validation-statistical-active-r20-v14)
                set_tests_properties(atcg3d_statistical_active_r20_validation PROPERTIES LABELS validation TIMEOUT 1800)
                add_test(NAME atcg3d_shared_angiogenesis_validation
                    COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/scripts/validate_shared_angiogenesis.py
                        --exe $<TARGET_FILE:atcg3d_shared_abm>
                        --config ${CMAKE_CURRENT_SOURCE_DIR}/ATCG3D_SharedRules/config/angiogenesis_controlled_v15.yaml
                        --seeds 16
                        --output ${CMAKE_CURRENT_BINARY_DIR}/validation-shared-angiogenesis)
                set_tests_properties(atcg3d_shared_angiogenesis_validation PROPERTIES LABELS validation TIMEOUT 1800)
                if(ATCG_ENABLE_LONG_VALIDATION)
                    foreach(long_case IN ITEMS r20_512_360 r200_512_360 r20_2000_720 r200_2000_720)
                        add_test(NAME atcg3d_long_${long_case}_validation
                            COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/scripts/validate_abm_vs_pde.py
                                --exe $<TARGET_FILE:atcg3d_shared_abm>
                                --config ${CMAKE_CURRENT_SOURCE_DIR}/ATCG3D_SharedRules/config/long_${long_case}_v14.yaml
                                --seeds 16 --sample-every-hours 24 --require-interior
                                --output ${CMAKE_CURRENT_BINARY_DIR}/validation-long-${long_case})
                        set_tests_properties(atcg3d_long_${long_case}_validation PROPERTIES LABELS validation_long TIMEOUT 604800)
                    endforeach()
                endif()
                add_test(NAME atcg3d_hybrid_invasion_validation
                    COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/scripts/validate_hybrid_invasion.py
                        --exe $<TARGET_FILE:atcg3d_hybrid>
                        --native-exe $<TARGET_FILE:atcg3d_shared_abm>
                        --seeds 16 --output ${CMAKE_CURRENT_BINARY_DIR}/validation-hybrid-invasion)
                set_tests_properties(atcg3d_hybrid_invasion_validation PROPERTIES LABELS validation TIMEOUT 3600)
                add_test(NAME atcg3d_hybrid_validation
                    COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/scripts/validate_hybrid.py
                        --exe $<TARGET_FILE:atcg3d_hybrid>
                        --config ${CMAKE_CURRENT_SOURCE_DIR}/ATCG3D_Hybrid/config/hybrid_regular_cycle_v2.yaml
                        --seeds 16 --output ${CMAKE_CURRENT_BINARY_DIR}/validation-hybrid)
                set_tests_properties(atcg3d_hybrid_validation PROPERTIES LABELS validation TIMEOUT 1800)
                add_test(NAME atcg3d_hybrid_volume_validation
                    COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/scripts/validate_hybrid.py
                        --exe $<TARGET_FILE:atcg3d_hybrid>
                        --config ${CMAKE_CURRENT_SOURCE_DIR}/ATCG3D_Hybrid/config/hybrid_regular_cycle_v3.yaml
                        --seeds 16 --output ${CMAKE_CURRENT_BINARY_DIR}/validation-hybrid-volume)
                set_tests_properties(atcg3d_hybrid_volume_validation PROPERTIES LABELS validation TIMEOUT 1800)
                add_test(NAME atcg3d_abm_pde_validation
                    COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/scripts/validate_abm_vs_pde.py
                        --exe $<TARGET_FILE:atcg3d_shared_abm>
                        --config ${CMAKE_CURRENT_SOURCE_DIR}/ATCG3D_SharedRules/config/validation_v7.yaml
                        --smoke-only --seeds 16 --output ${CMAKE_CURRENT_BINARY_DIR}/validation)
                set_tests_properties(atcg3d_abm_pde_validation PROPERTIES LABELS validation TIMEOUT 180)
                add_test(NAME atcg3d_regular_cycle_validation
                    COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/scripts/validate_abm_vs_pde.py
                        --exe $<TARGET_FILE:atcg3d_shared_abm>
                        --config ${CMAKE_CURRENT_SOURCE_DIR}/ATCG3D_SharedRules/config/regular_cycle_v10.yaml
                        --smoke-only --seeds 16 --output ${CMAKE_CURRENT_BINARY_DIR}/validation-regular)
                set_tests_properties(atcg3d_regular_cycle_validation PROPERTIES LABELS validation TIMEOUT 300)
                add_test(NAME atcg3d_vascular_validation
                    COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/scripts/validate_abm_vs_pde.py
                        --exe $<TARGET_FILE:atcg3d_shared_abm>
                        --config ${CMAKE_CURRENT_SOURCE_DIR}/ATCG3D_SharedRules/config/angiogenesis_v8.yaml
                        --smoke-only --seeds 16 --vascular --output ${CMAKE_CURRENT_BINARY_DIR}/validation-vascular)
                set_tests_properties(atcg3d_vascular_validation PROPERTIES LABELS validation TIMEOUT 180)
                foreach(multiplier IN ITEMS 20 200)
                    if(multiplier EQUAL 200)
                        set(active_config active_r200_ci_v12.yaml)
                    else()
                        set(active_config active_r20_v12.yaml)
                    endif()
                    add_test(NAME atcg3d_active_r${multiplier}_validation
                        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/scripts/validate_abm_vs_pde.py
                            --exe $<TARGET_FILE:atcg3d_shared_abm>
                            --config ${CMAKE_CURRENT_SOURCE_DIR}/ATCG3D_SharedRules/config/${active_config}
                            --smoke-only --seeds 16 --minimum-active-fraction 0.1
                            --output ${CMAKE_CURRENT_BINARY_DIR}/validation-active-r${multiplier})
                    set_tests_properties(atcg3d_active_r${multiplier}_validation PROPERTIES LABELS validation TIMEOUT 1800)
                endforeach()
                add_test(NAME atcg3d_nonlinear_vascular_validation
                    COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/scripts/validate_abm_vs_pde.py
                        --exe $<TARGET_FILE:atcg3d_shared_abm>
                        --config ${CMAKE_CURRENT_SOURCE_DIR}/ATCG3D_SharedRules/config/angiogenesis_nonlinear_v8.yaml
                        --smoke-only --seeds 16 --vascular --require-vascular-activity
                        --output ${CMAKE_CURRENT_BINARY_DIR}/validation-vascular-nonlinear)
                set_tests_properties(atcg3d_nonlinear_vascular_validation PROPERTIES LABELS validation TIMEOUT 180)
            endif()
            add_executable(atcg3d_structured_pde_test
                tests/3d/structured_pde_test.cpp)
            target_link_libraries(atcg3d_structured_pde_test PRIVATE
                atcg3d_structured_pde_core)
            target_compile_options(atcg3d_structured_pde_test PRIVATE -UNDEBUG)
            target_compile_definitions(atcg3d_structured_pde_test PRIVATE
                ATCG_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")
            add_test(NAME atcg3d_structured_pde_test
                COMMAND atcg3d_structured_pde_test)
        endif()
        if(ATCG3D_ENABLE_HDF5_CHECKPOINT)
            add_executable(atcg3d_checkpoint_test tests/3d/checkpoint_test.cpp)
            target_link_libraries(atcg3d_checkpoint_test PRIVATE atcg3d_core)
            target_compile_options(atcg3d_checkpoint_test PRIVATE -UNDEBUG)
            add_test(NAME atcg3d_checkpoint_test COMMAND atcg3d_checkpoint_test)
        endif()
        if(ATCG_BUILD_3D_NUTRIENT AND ATCG_BUILD_3D_STRUCTURED_PDE)
            add_executable(atcg3d_compatibility_test tests/3d/compatibility_test.cpp)
            target_link_libraries(atcg3d_compatibility_test PRIVATE
                atcg3d_nutrient_core atcg3d_structured_pde_core ${ATCG3D_YAML_TARGET})
            target_compile_options(atcg3d_compatibility_test PRIVATE -UNDEBUG)
            target_compile_definitions(atcg3d_compatibility_test PRIVATE
                ATCG_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")
            add_test(NAME atcg3d_compatibility_test COMMAND atcg3d_compatibility_test)
            foreach(contract_test IN ITEMS alignment legacy_pde_regression)
                add_executable(atcg3d_${contract_test}_test tests/3d/${contract_test}_test.cpp)
                target_link_libraries(atcg3d_${contract_test}_test PRIVATE
                    atcg3d_nutrient_core atcg3d_structured_pde_core)
                target_compile_options(atcg3d_${contract_test}_test PRIVATE -UNDEBUG)
                target_compile_definitions(atcg3d_${contract_test}_test PRIVATE
                    ATCG_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")
                add_test(NAME atcg3d_${contract_test}_test COMMAND atcg3d_${contract_test}_test)
            endforeach()
        endif()
        if(ATCG3D_ENABLE_VTKHDF)
            add_executable(atcg3d_vtkhdf_test tests/3d/vtkhdf_test.cpp)
            target_link_libraries(atcg3d_vtkhdf_test PRIVATE atcg3d_core)
            target_compile_options(atcg3d_vtkhdf_test PRIVATE -UNDEBUG)
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
        add_test(
            NAME atcg3d_angiogenesis_scale_smoke
            COMMAND ${CMAKE_COMMAND}
                -DBENCHMARK=$<TARGET_FILE:atcg3d_angiogenesis_scale_benchmark>
                -DCELLS=1000
                -DEVENTS=64
                -P ${CMAKE_CURRENT_SOURCE_DIR}/tests/3d/angiogenesis_scale_test.cmake
        )
        set_tests_properties(atcg3d_angiogenesis_scale_smoke PROPERTIES
            LABELS benchmark
            PROCESSORS 4
            TIMEOUT 90
        )
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
            set_tests_properties(atcg3d_ten_million_core_benchmark PROPERTIES
                LABELS large
                RUN_SERIAL TRUE
            )
            add_test(
                NAME atcg3d_ten_million_angiogenesis_benchmark
                COMMAND atcg3d_angiogenesis_scale_benchmark
                    --cells 10000000
                    --events 256
                    --threads 1
            )
            set_tests_properties(
                atcg3d_ten_million_angiogenesis_benchmark
                PROPERTIES
                    LABELS large
                    RUN_SERIAL TRUE
            )
            add_test(
                NAME atcg3d_ten_million_angiogenesis_production_benchmark
                COMMAND ${CMAKE_COMMAND}
                    -DBENCHMARK=$<TARGET_FILE:atcg3d_angiogenesis_scale_benchmark>
                    -DCELLS=10000000
                    -DEVENTS=100000
                    -DPROFILE=production
                    -P ${CMAKE_CURRENT_SOURCE_DIR}/tests/3d/angiogenesis_scale_test.cmake
            )
            set_tests_properties(
                atcg3d_ten_million_angiogenesis_production_benchmark
                PROPERTIES
                    LABELS large
                    RUN_SERIAL TRUE
                    PROCESSORS 4
                    TIMEOUT 1800
            )
        endif()
    endif()
endif()
