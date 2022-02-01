# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.1...3.22)
project(Orthtree_Examples)

find_package(CGAL REQUIRED QUIET OPTIONAL_COMPONENTS Core)
if (CGAL_FOUND)

    create_single_source_cgal_program("octree_build_from_point_set.cpp")
    create_single_source_cgal_program("octree_build_from_point_vector.cpp")
    create_single_source_cgal_program("octree_build_with_custom_split.cpp")
    create_single_source_cgal_program("octree_find_nearest_neighbor.cpp")
    create_single_source_cgal_program("octree_traversal_custom.cpp")
    create_single_source_cgal_program("octree_traversal_manual.cpp")
    create_single_source_cgal_program("octree_traversal_preorder.cpp")
    create_single_source_cgal_program("octree_grade.cpp")
    create_single_source_cgal_program("quadtree_build_from_point_vector.cpp")

    find_package(Eigen3 3.1.91) #(requires 3.2.0 or greater)
    include(CGAL_Eigen_support)
    if (TARGET CGAL::Eigen_support)
      create_single_source_cgal_program("orthtree_build.cpp")
      target_link_libraries(orthtree_build PUBLIC CGAL::Eigen_support)
    endif()

else ()
    message(WARNING
            "This program requires the CGAL library, and will not be compiled.")
endif ()
