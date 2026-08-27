cmake_minimum_required(VERSION 3.12...3.31)
project(Mesh_smoothing_3_Examples)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

message(STATUS "CGAL_DIR: ${CGAL_DIR}")

find_package(CGAL REQUIRED)
find_package(Eigen3 QUIET)

if (NOT TARGET CGAL::Eigen3_support)
  message("NOTICE: This project requires the Eigen library, and will not be compiled.")
  return()
endif()

create_single_source_cgal_program("c3t3_smooth.cpp")
target_link_libraries(c3t3_smooth PRIVATE CGAL::Eigen3_support)

find_package(OpenMP COMPONENTS CXX)
if(OpenMP_CXX_FOUND)
    message(STATUS "OpenMP found")
    if(MSVC)
        # Override the compile options for the imported OpenMP target
        set_target_properties(OpenMP::OpenMP_CXX PROPERTIES
                            INTERFACE_COMPILE_OPTIONS "/openmp:llvm")
    endif()
    # target_link_libraries(alpha_wrap_improvement PRIVATE OpenMP::OpenMP_CXX)
    target_link_libraries(c3t3_smooth PRIVATE OpenMP::OpenMP_CXX)
else()
    message(STATUS "OpenMP was not found.")
endif()


option(MAKE_LIB_LINKING_EXAMPLES "Build examples of other meshing libs that link against Mesh_smoothing_3" OFF)

if (MAKE_LIB_LINKING_EXAMPLES)
    create_single_source_cgal_program("alpha_wrap_improvement.cpp")
    target_link_libraries(alpha_wrap_improvement PRIVATE CGAL::Eigen3_support)
    if(OpenMP_CXX_FOUND)
        target_link_libraries(alpha_wrap_improvement PRIVATE OpenMP::OpenMP_CXX)
    endif()

    include(FetchContent)
    add_subdirectory(Cinolib)
    add_subdirectory(libigl)
    add_subdirectory(geogram)
endif()
