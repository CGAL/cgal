cmake_minimum_required(VERSION 3.12...3.31)
project(Kernel_23_Tests)

# CGAL

find_package(CGAL REQUIRED COMPONENTS Core)
IF ( CGAL_FOUND )
	message("${CGAL_LIBRARY_DIRS}")
ENDIF()

include_directories(BEFORE "include")

# For Tetgen
#set(CMAKE_TETGEN_PATH "/Users/umenohana/Dropbox/G-Mod/TopAlg/code/tetgen1.6.0/build/" )
#configure_file( ${CGAL_LIBRARY_DIRS}/include/CGAL/HDVF/config.h.in ${PROJECT_SOURCE_DIR}/generated/config.h )
#include_directories( ${PROJECT_SOURCE_DIR}/include/CGAL/HDVF/generated/ )

add_executable(hdvf main_hdvf.cpp arguments.cpp)
target_link_libraries(hdvf CGAL::CGAL CGAL::CGAL_Core)

add_executable(per_hdvf main_per_hdvf.cpp arguments.cpp)
target_link_libraries(per_hdvf CGAL::CGAL CGAL::CGAL_Core)

add_executable(dual_hdvf main_dual_hdvf.cpp arguments.cpp)
target_link_libraries(dual_hdvf CGAL::CGAL CGAL::CGAL_Core)

add_executable(test_hdvf test_hdvf_example.cpp)
target_link_libraries(test_hdvf CGAL::CGAL CGAL::CGAL_Core)

add_executable(example_hdvf_simplicial example_hdvf_simplicial.cpp)
target_link_libraries(example_hdvf_simplicial CGAL::CGAL CGAL::CGAL_Core)

add_executable(example_hdvf_cubical example_hdvf_cubical.cpp)
target_link_libraries(example_hdvf_cubical CGAL::CGAL CGAL::CGAL_Core)

add_executable(example_matrix_chain example_matrix_chain.cpp)
target_link_libraries(example_matrix_chain CGAL::CGAL CGAL::CGAL_Core)
