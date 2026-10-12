cmake_minimum_required(VERSION 3.12...3.31)
project(Feature_graph_Examples)

add_compile_definitions(CGAL_MESH_3_NO_DEPRECATED_SURFACE_INDEX
                        CGAL_MESH_3_NO_DEPRECATED_C3T3_ITERATORS)

if(CGAL_MESH_3_VERBOSE)
  add_compile_definitions(CGAL_MESH_3_VERBOSE)
endif()

find_package(CGAL REQUIRED COMPONENTS ImageIO)

option(CGAL_ACTIVATE_CONCURRENT_MESH_3 "Activate parallelism in Mesh_3" OFF)

if(CGAL_ACTIVATE_CONCURRENT_MESH_3 OR "$ENV{CGAL_ACTIVATE_CONCURRENT_MESH_3}")
  add_compile_definitions(CGAL_CONCURRENT_MESH_3)
  find_package(TBB REQUIRED)
  include(CGAL_TBB_support)
endif()

# Use Eigen
find_package(Eigen3 QUIET)
include(CGAL_Eigen3_support)

if(NOT TARGET CGAL::Eigen3_support)
  message("NOTICE: All examples require the Eigen3 library, and will not be compiled.")
  return()
endif()

# Use DGtal
find_package(DGtal QUIET)
if(DGtal_FOUND AND NOT TARGET CGAL::DGTAL_support)
  add_library(CGAL::DGTAL_support INTERFACE IMPORTED)
  set_target_properties(CGAL::DGTAL_support PROPERTIES
  	INTERFACE_COMPILE_DEFINITIONS "CGAL_DGTAL_ENABLED"
    INTERFACE_INCLUDE_DIRECTORIES "${DGTAL_INCLUDE_DIRS}"
    INTERFACE_LINK_LIBRARIES "${DGTAL_LIBRARIES}")
  set(CGAL_Feature_graph_USE_DGTAL ON CACHE BOOL "CGAL_Feature_graph_USE_DGTAL uses DGtal")
endif()

if(TARGET CGAL::CGAL_ImageIO)

  if(CGAL_ImageIO_USE_ZLIB)

    if(CGAL_Feature_graph_USE_DGTAL)
      create_single_source_cgal_program(
        "feature_graph_from_segmented_image.cpp")
      target_link_libraries(feature_graph_from_segmented_image
                            PUBLIC CGAL::Eigen3_support CGAL::DGTAL_support)
    else()
      message(STATUS "NOTICE: The example feature_graph_from_segmented_image.cpp need CGAL_Feature_graph to be configured with DGtal support, and will not be compiled.")
    endif()
  else()
    message(STATUS "NOTICE: The example feature_graph_from_segmented_image.cpp need CGAL_ImageIO to be configured with ZLIB support, and will not be compiled.")
  endif()

else()
  message(STATUS "NOTICE: Some examples need the CGAL_ImageIO library, and will not be compiled.")
endif()

if(CGAL_ACTIVATE_CONCURRENT_MESH_3 AND TARGET CGAL::TBB_support)
  foreach(
    target
    feature_graph_from_segmented_image)
    if(TARGET ${target})
      target_link_libraries(${target} PRIVATE CGAL::TBB_support)
    endif()
  endforeach()
endif()
