if(NOT DEFINED BENCHMARK OR BENCHMARK STREQUAL "")
    message(FATAL_ERROR "BENCHMARK executable was not provided")
endif()
if(NOT DEFINED CELLS)
    set(CELLS 1000)
endif()
if(NOT DEFINED EVENTS)
    set(EVENTS 64)
endif()
if(NOT DEFINED PROFILE)
    set(PROFILE synthetic)
endif()
if(NOT PROFILE STREQUAL "synthetic" AND NOT PROFILE STREQUAL "production")
    message(FATAL_ERROR "PROFILE must be synthetic or production")
endif()

function(run_case thread_count output_name)
    execute_process(
        COMMAND "${BENCHMARK}"
            --cells "${CELLS}"
            --events "${EVENTS}"
            --threads "${thread_count}"
            --profile "${PROFILE}"
        RESULT_VARIABLE result
        OUTPUT_VARIABLE output
        ERROR_VARIABLE error
        TIMEOUT 60
    )
    if(NOT result STREQUAL "0")
        message(FATAL_ERROR
            "angiogenesis benchmark failed for ${thread_count} threads "
            "(exit ${result})\nstdout:\n${output}\nstderr:\n${error}")
    endif()
    set(${output_name} "${output}" PARENT_SCOPE)
endfunction()

function(require_exact output field expected)
    string(REGEX MATCH
        "\"${field}\"[ \t\r\n]*:[ \t\r\n]*([0-9][0-9]*)"
        matched "${output}")
    if(matched STREQUAL "" OR NOT CMAKE_MATCH_1 STREQUAL "${expected}")
        message(FATAL_ERROR
            "expected ${field}=${expected}; benchmark output was:\n${output}")
    endif()
endfunction()

function(require_positive output field)
    string(REGEX MATCH
        "\"${field}\"[ \t\r\n]*:[ \t\r\n]*([1-9][0-9]*)"
        matched "${output}")
    if(matched STREQUAL "")
        message(FATAL_ERROR
            "expected positive ${field}; benchmark output was:\n${output}")
    endif()
endfunction()

function(extract_checksum output result_name)
    string(REGEX MATCH
        "\"checksum\"[ \t\r\n]*:[ \t\r\n]*([0-9][0-9]*)"
        matched "${output}")
    if(matched STREQUAL "")
        message(FATAL_ERROR "benchmark output has no checksum:\n${output}")
    endif()
    set(${result_name} "${CMAKE_MATCH_1}" PARENT_SCOPE)
endfunction()

run_case(1 one_thread)
run_case(4 four_threads)

foreach(output IN ITEMS one_thread four_threads)
    if(PROFILE STREQUAL "synthetic")
        require_exact("${${output}}" seed_attempts 1)
        require_exact("${${output}}" committed_roots 1)
        require_exact("${${output}}" seed_rejections 0)
        require_exact("${${output}}" root_nodes 1)
        require_exact("${${output}}" vessel_tips 2)
    else()
        require_positive("${${output}}" seed_attempts)
        require_positive("${${output}}" committed_roots)
        require_positive("${${output}}" root_nodes)
        require_positive("${${output}}" vessel_tips)
    endif()
    foreach(field IN ITEMS
            vessel_growth_attempts vessel_growth_commits
            vascular_displacements inward_nodes outward_nodes
            vessel_occupied_voxels vascular_influenced_voxels
            tracked_component_bytes peak_rss_bytes)
        require_positive("${${output}}" "${field}")
    endforeach()
endforeach()

extract_checksum("${one_thread}" one_checksum)
extract_checksum("${four_threads}" four_checksum)
if(NOT one_checksum STREQUAL four_checksum)
    message(FATAL_ERROR
        "1/4-thread angiogenesis checksums differ: "
        "${one_checksum} != ${four_checksum}\n"
        "one thread:\n${one_thread}\nfour threads:\n${four_threads}")
endif()

message(STATUS
    "angiogenesis scale passed: profile=${PROFILE}, cells=${CELLS}, events=${EVENTS}, "
    "checksum=${one_checksum}")
