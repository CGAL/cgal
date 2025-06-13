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
# Fix the metadata generation command
# Path to your Python script - note the path correction
set(METADATA_SCRIPT "${CMAKE_CURRENT_SOURCE_DIR}/generate_source_exec_metadata.py")

#add_compile_definitions(CGAL_MESH_3_NO_DEPRECATED_SURFACE_INDEX
#               -DCGAL_MESH_3_NO_DEPRECATED_C3T3_ITERATORS)

# Activate verbose mode? (turned OFF by default)
option(CGAL_ACTIVATE_CONCURRENT_MESH_3 "Activate verbose mode in Mesh_3" OFF)
if(MESH_3_VERBOSE)
  add_compile_definitions(CGAL_MESH_3_VERBOSE)
endif()

# Activate concurrency? (turned OFF by default)
option(CGAL_ACTIVATE_CONCURRENT_MESH_3 "Activate parallelism in Mesh_3" OFF)

# And add -DCGAL_CONCURRENT_MESH_3 if that option is ON
if(CGAL_ACTIVATE_CONCURRENT_MESH_3)
  add_compile_definitions(CGAL_CONCURRENT_MESH_3)
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
create_single_source_cgal_program("benchmark_atomic_tetrahedral_remeshing.cpp")

target_include_directories(benchmark_tetrahedral_remeshing PRIVATE ${CMAKE_SOURCE_DIR}/Benchmark_creation/external/)
target_include_directories(benchmark_adaptive_tetrahedral_remeshing PRIVATE ${CMAKE_SOURCE_DIR}/Benchmark_creation/external)
target_include_directories(benchmark_atomic_tetrahedral_remeshing PRIVATE ${CMAKE_SOURCE_DIR}/Benchmark_creation/external)

if(TARGET CGAL::TBB_support)
  target_link_libraries(benchmark_tetrahedral_remeshing PRIVATE CGAL::TBB_support)
  target_link_libraries(benchmark_adaptive_tetrahedral_remeshing PRIVATE CGAL::TBB_support)
  target_link_libraries(benchmark_atomic_tetrahedral_remeshing PRIVATE CGAL::TBB_support)
endif()

# Link with Boost.ProgramOptions (optional)
find_package(Boost QUIET COMPONENTS program_options)
if(Boost_PROGRAM_OPTIONS_FOUND)
  target_link_libraries(benchmark_tetrahedral_remeshing PRIVATE Boost::program_options)
  target_link_libraries(benchmark_adaptive_tetrahedral_remeshing PRIVATE Boost::program_options)
  target_link_libraries(benchmark_atomic_tetrahedral_remeshing PRIVATE Boost::program_options)
endif()

# Add a post-build step to generate the metadata file
add_custom_command(TARGET benchmark_tetrahedral_remeshing POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E echo "Generating metadata for $<TARGET_FILE:benchmark_tetrahedral_remeshing>"
    # Combine environment setting and Python execution into a single command
    COMMAND ${CMAKE_COMMAND} -E env 
            "PYTHONPATH=${CMAKE_SOURCE_DIR}"
            "${Python_EXECUTABLE}" 
            "${METADATA_SCRIPT}"
            "$<TARGET_FILE:benchmark_tetrahedral_remeshing>"
            "benchmark_tetrahedral_remeshing.cpp"
    WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
    COMMENT "Generating metadata file for benchmark_tetrahedral_remeshing"
    VERBATIM
)

add_custom_command(TARGET benchmark_adaptive_tetrahedral_remeshing POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E echo "Generating metadata for $<TARGET_FILE:benchmark_adaptive_tetrahedral_remeshing>"
    COMMAND ${CMAKE_COMMAND} -E env 
            "PYTHONPATH=${CMAKE_SOURCE_DIR}"
            "${Python_EXECUTABLE}" 
            "${METADATA_SCRIPT}"
            "$<TARGET_FILE:benchmark_adaptive_tetrahedral_remeshing>"
            "benchmark_adaptive_tetrahedral_remeshing.cpp"
    WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
    COMMENT "Generating metadata file for benchmark_adaptive_tetrahedral_remeshing"
    VERBATIM
)

add_custom_command(TARGET benchmark_atomic_tetrahedral_remeshing POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E echo "Generating metadata for $<TARGET_FILE:benchmark_atomic_tetrahedral_remeshing>"
    COMMAND ${CMAKE_COMMAND} -E env 
            "PYTHONPATH=${CMAKE_SOURCE_DIR}"
            "${Python_EXECUTABLE}" 
            "${METADATA_SCRIPT}"
            "$<TARGET_FILE:benchmark_atomic_tetrahedral_remeshing>"
            "benchmark_atomic_tetrahedral_remeshing.cpp"
    WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
    COMMENT "Generating metadata file for benchmark_atomic_tetrahedral_remeshing"
    VERBATIM
)

