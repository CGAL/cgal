# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.12...3.31)
project(Tetrahedral_remeshing_benchmark)

find_package(CGAL REQUIRED)

find_package(Python 3 COMPONENTS Interpreter REQUIRED)
if (Python_FOUND)
    message(STATUS "Found Python interpreter: ${Python_EXECUTABLE}")
else()
    message(FATAL_ERROR "Python interpreter not found. Please install Python 3 or ensure it's in your PATH.")
endif()
# Prefer a usable Python executable (Windows Store paths may be non-invokable). Build a fallback.
set(METADATA_PYTHON "")
if (EXISTS "${Python_EXECUTABLE}")
  # Usable absolute path
  set(METADATA_PYTHON "${Python_EXECUTABLE}")
else()
  # Try the Windows 'py' launcher first
  find_program(PY_LAUNCHER py)
  if (PY_LAUNCHER)
    set(METADATA_PYTHON "${PY_LAUNCHER}" -3)
  else()
    # Fallback to 'python' on PATH
    find_program(SYS_PYTHON python)
    if (SYS_PYTHON)
      set(METADATA_PYTHON "${SYS_PYTHON}")
    else()
      message(FATAL_ERROR "No usable Python interpreter found (Python_EXECUTABLE invalid, no 'py' launcher, and no 'python' on PATH).")
    endif()
  endif()
endif()
# Fix the metadata generation command
# Path to your Python script - note the path correction
set(METADATA_SCRIPT "${CMAKE_CURRENT_SOURCE_DIR}/generate_source_exec_metadata.py")

option(CGAL_ACTIVATE_CONCURRENT_REMESHING "Activate parallelism in Tetrahedral_remeshig" ON)

# Enable debug output for edge flip operations only in debug mode
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
  add_compile_definitions(CGAL_REFACTORED_TETRAHEDRAL_REMESHING_DEBUG)
endif()

if(CGAL_ACTIVATE_CONCURRENT_REMESHING)
  message(STATUS "Activating concurrent remeshing")
  #add_compile_definitions(CGAL_CONCURRENT_TETRAHEDRAL_REMESHING)
  find_package(TBB REQUIRED)
  include(CGAL_TBB_support)
else()
  option(LINK_WITH_TBB "Link with TBB anyway so we can use TBB timers for profiling" ON)
  if(LINK_WITH_TBB)
    find_package(TBB)
    include(CGAL_TBB_support)
  endif(LINK_WITH_TBB)
endif()

# Compilable benchmark
create_single_source_cgal_program("benchmark_tetrahedral_remeshing.cpp")
create_single_source_cgal_program("benchmark_adaptive_tetrahedral_remeshing.cpp")
create_single_source_cgal_program("benchmark_refactored_tetrahedral_remeshing.cpp")
create_single_source_cgal_program("benchmark_original_tetrahedral_remeshing.cpp")
create_single_source_cgal_program("benchmark_refactored_sequential_tetrahedral_remeshing.cpp")
create_single_source_cgal_program("benchmark_domain_decomp_extract_merge_tetrahedral_remeshing.cpp")
create_single_source_cgal_program("make_corpus.cpp")
create_single_source_cgal_program("benchmark_domain_decomp_seq_global_pass.cpp")

target_include_directories(benchmark_tetrahedral_remeshing PRIVATE ${CMAKE_SOURCE_DIR}/Benchmark_creation/external/)
target_include_directories(benchmark_adaptive_tetrahedral_remeshing PRIVATE ${CMAKE_SOURCE_DIR}/Benchmark_creation/external)
target_include_directories(benchmark_refactored_tetrahedral_remeshing PRIVATE ${CMAKE_SOURCE_DIR}/Benchmark_creation/external)
target_include_directories(benchmark_original_tetrahedral_remeshing PRIVATE ${CMAKE_SOURCE_DIR}/Benchmark_creation/external)
target_include_directories(benchmark_refactored_sequential_tetrahedral_remeshing PRIVATE ${CMAKE_SOURCE_DIR}/Benchmark_creation/external)
target_include_directories(benchmark_domain_decomp_extract_merge_tetrahedral_remeshing PRIVATE ${CMAKE_SOURCE_DIR}/Benchmark_creation/external)

# METIS for the ParMmg-style dual-graph partitioner in the domain-decomposition
# benchmark. Optional: without it the benchmark falls back to the geometric
# grid partition (CGAL_DD_USE_METIS undefined).
find_path(METIS_INCLUDE_DIR metis.h)
find_library(METIS_LIBRARY NAMES metis)
find_library(GKLIB_LIBRARY NAMES GKlib gklib)
if(METIS_INCLUDE_DIR AND METIS_LIBRARY)
  message(STATUS "DD benchmark: METIS found at ${METIS_LIBRARY}")
  target_include_directories(benchmark_domain_decomp_extract_merge_tetrahedral_remeshing PRIVATE ${METIS_INCLUDE_DIR})
  target_compile_definitions(benchmark_domain_decomp_extract_merge_tetrahedral_remeshing PRIVATE CGAL_DD_USE_METIS)
  target_link_libraries(benchmark_domain_decomp_extract_merge_tetrahedral_remeshing PRIVATE ${METIS_LIBRARY})
  if(GKLIB_LIBRARY)
    target_link_libraries(benchmark_domain_decomp_extract_merge_tetrahedral_remeshing PRIVATE ${GKLIB_LIBRARY})
  endif()
else()
  message(STATUS "DD benchmark: METIS not found, using the geometric grid partition")
endif()

# ParMETIS for the adaptive REpartitioner (ParMETIS_V3_AdaptiveRepart). Also
# optional: without it the benchmark falls back to from-scratch serial METIS.
# It is driven with a single MPI rank -- the parallelism here is TBB -- but it
# is still an MPI library, so MPI must be linked and MPI_Init called.
find_path(PARMETIS_INCLUDE_DIR parmetis.h)
find_library(PARMETIS_LIBRARY NAMES parmetis)
find_package(MPI QUIET COMPONENTS C)
if(PARMETIS_INCLUDE_DIR AND PARMETIS_LIBRARY AND MPI_C_FOUND)
  message(STATUS "DD benchmark: ParMETIS found at ${PARMETIS_LIBRARY}")
  target_include_directories(benchmark_domain_decomp_extract_merge_tetrahedral_remeshing PRIVATE ${PARMETIS_INCLUDE_DIR})
  target_compile_definitions(benchmark_domain_decomp_extract_merge_tetrahedral_remeshing PRIVATE CGAL_DD_USE_PARMETIS)
  target_link_libraries(benchmark_domain_decomp_extract_merge_tetrahedral_remeshing PRIVATE ${PARMETIS_LIBRARY} MPI::MPI_C)
else()
  message(STATUS "DD benchmark: ParMETIS not found, adaptive repartitioning disabled")
endif()
target_include_directories(benchmark_domain_decomp_seq_global_pass PRIVATE ${CMAKE_SOURCE_DIR}/Benchmark_creation/external)

if(TARGET CGAL::TBB_support)
  target_link_libraries(benchmark_tetrahedral_remeshing PRIVATE CGAL::TBB_support)
  target_link_libraries(benchmark_adaptive_tetrahedral_remeshing PRIVATE CGAL::TBB_support)
  target_link_libraries(benchmark_refactored_tetrahedral_remeshing PRIVATE CGAL::TBB_support)
  target_link_libraries(benchmark_original_tetrahedral_remeshing PRIVATE CGAL::TBB_support)
  target_link_libraries(benchmark_refactored_sequential_tetrahedral_remeshing PRIVATE CGAL::TBB_support)
  target_link_libraries(benchmark_domain_decomp_extract_merge_tetrahedral_remeshing PRIVATE CGAL::TBB_support)
  target_link_libraries(benchmark_domain_decomp_seq_global_pass PRIVATE CGAL::TBB_support)
endif()

# Link with Boost.ProgramOptions (optional)
find_package(Boost QUIET COMPONENTS program_options)
if(Boost_PROGRAM_OPTIONS_FOUND)
  target_link_libraries(benchmark_tetrahedral_remeshing PRIVATE Boost::program_options)
  target_link_libraries(benchmark_adaptive_tetrahedral_remeshing PRIVATE Boost::program_options)
  target_link_libraries(benchmark_refactored_tetrahedral_remeshing PRIVATE Boost::program_options)
  target_link_libraries(benchmark_original_tetrahedral_remeshing PRIVATE Boost::program_options)
  target_link_libraries(benchmark_refactored_sequential_tetrahedral_remeshing PRIVATE Boost::program_options)
  target_link_libraries(benchmark_domain_decomp_extract_merge_tetrahedral_remeshing PRIVATE Boost::program_options)
  target_link_libraries(benchmark_domain_decomp_seq_global_pass PRIVATE Boost::program_options)
endif()

# Add a post-build step to generate the metadata file
add_custom_command(TARGET benchmark_tetrahedral_remeshing POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E echo "Generating metadata for $<TARGET_FILE:benchmark_tetrahedral_remeshing>"
    # Combine environment setting and Python execution into a single command
    COMMAND ${CMAKE_COMMAND} -E env
            "PYTHONPATH=${CMAKE_CURRENT_SOURCE_DIR}"
            ${METADATA_PYTHON}
            "${METADATA_SCRIPT}"
            "$<TARGET_FILE:benchmark_tetrahedral_remeshing>"
            "${CMAKE_CURRENT_SOURCE_DIR}/benchmark_tetrahedral_remeshing.cpp"
    WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
    COMMENT "Generating metadata file for benchmark_tetrahedral_remeshing"
    VERBATIM
)

add_custom_command(TARGET benchmark_adaptive_tetrahedral_remeshing POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E echo "Generating metadata for $<TARGET_FILE:benchmark_adaptive_tetrahedral_remeshing>"
    COMMAND ${CMAKE_COMMAND} -E env
            "PYTHONPATH=${CMAKE_CURRENT_SOURCE_DIR}"
            ${METADATA_PYTHON}
            "${METADATA_SCRIPT}"
            "$<TARGET_FILE:benchmark_adaptive_tetrahedral_remeshing>"
            "${CMAKE_CURRENT_SOURCE_DIR}/benchmark_adaptive_tetrahedral_remeshing.cpp"
    WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
    COMMENT "Generating metadata file for benchmark_adaptive_tetrahedral_remeshing"
    VERBATIM
)

add_custom_command(TARGET benchmark_refactored_tetrahedral_remeshing POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E echo "Generating metadata for $<TARGET_FILE:benchmark_refactored_tetrahedral_remeshing>"
    COMMAND ${CMAKE_COMMAND} -E env
            "PYTHONPATH=${CMAKE_CURRENT_SOURCE_DIR}"
            ${METADATA_PYTHON}
            "${METADATA_SCRIPT}"
            "$<TARGET_FILE:benchmark_refactored_tetrahedral_remeshing>"
            "${CMAKE_CURRENT_SOURCE_DIR}/benchmark_refactored_tetrahedral_remeshing.cpp"
    WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
    COMMENT "Generating metadata file for benchmark_refactored_tetrahedral_remeshing"
    VERBATIM
)

