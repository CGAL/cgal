# Created by the script cgal_create_CMakeLists
# This is the CMake script for compiling a set of CGAL applications.

cmake_minimum_required(VERSION 3.12...3.31)
project(Classification_Examples)

# CGAL and its components
find_package(CGAL REQUIRED)

set(Classification_dependencies_met TRUE)

find_package(Boost OPTIONAL_COMPONENTS serialization iostreams)
include(CGAL_Boost_serialization_support)
include(CGAL_Boost_iostreams_support)

if(NOT TARGET CGAL::Boost_serialization_support)
  message("NOTICE: This project requires Boost Serialization, and will not be compiled.")
  set(Classification_dependencies_met FALSE)
endif()

if(NOT TARGET CGAL::Boost_iostreams_support)
  message("NOTICE: This project requires Boost IO Streams, and will not be compiled.")
  set(Classification_dependencies_met FALSE)
endif()

find_package(Eigen3 3.1.0 QUIET) #(3.1.0 or greater)
include(CGAL_Eigen3_support)
if(NOT TARGET CGAL::Eigen3_support)
  message("NOTICE: This project requires the Eigen library, and will not be compiled.")
  set(Classification_dependencies_met FALSE)
endif()

if(NOT Classification_dependencies_met)
  return()
endif()

find_package(TBB QUIET)
include(CGAL_TBB_support)

create_single_source_cgal_program( "example_classification.cpp" )
create_single_source_cgal_program( "example_ethz_random_forest.cpp" )
create_single_source_cgal_program( "example_feature.cpp" )
create_single_source_cgal_program( "example_advanced_features.cpp" )
create_single_source_cgal_program( "example_generate_advanced_features.cpp" )
create_single_source_cgal_program( "example_multi_dim_feature.cpp" )
create_single_source_cgal_program( "example_generation_and_training.cpp" )
create_single_source_cgal_program( "example_mesh_classification.cpp" )
create_single_source_cgal_program( "example_cluster_classification.cpp" )
create_single_source_cgal_program( "gis_tutorial_example.cpp" )
create_single_source_cgal_program( "example_deprecated_conversion.cpp" )

find_package(OpenCV QUIET COMPONENTS core ml) # Need core + machine learning
include(CGAL_OpenCV_support)
if(TARGET CGAL::OpenCV_support)
  message(STATUS "Found OpenCV")

  create_single_source_cgal_program( "example_opencv_random_forest.cpp" )
  target_link_libraries(example_opencv_random_forest PRIVATE CGAL::OpenCV_support)
else()
  message("NOTICE: OpenCV was not found. OpenCV random forest predicate for classification won't be available.")
endif()

foreach(target
    example_classification
    example_ethz_random_forest
    example_feature
    example_advanced_features
    example_generate_advanced_features
    example_multi_dim_feature
    example_generation_and_training
    example_mesh_classification
    example_cluster_classification
    example_opencv_random_forest
    gis_tutorial_example
    example_deprecated_conversion)
  if(TARGET ${target})
    target_link_libraries(${target} PRIVATE CGAL::Eigen3_support
                                            CGAL::Boost_iostreams_support
                                            CGAL::Boost_serialization_support)
    if(TARGET CGAL::TBB_support)
      target_link_libraries(${target} PRIVATE CGAL::TBB_support)
    endif()
  endif()
endforeach()

# TODO: Add a proper CGAL::Embree_support *.cmake file
find_package(embree QUIET)
message("Embree found: ${embree_FOUND}")
include(CGAL_embree_support)
if(TARGET CGAL::Embree_support)#if (embree_FOUND)#if (EMBREE_VERSION)
message("Adding embree support...")
#  include_directories(${EMBREE_INCLUDE_DIRS})
#  include_directories(${EMBREE_DIR})
#  target_include_directories(example_advanced_features PUBLIC ${EMBREE_INCLUDE_DIRS})
#  target_link_libraries(example_advanced_features PUBLIC ${EMBREE_LIBRARIES})
#target_compile_definitions(example_advanced_features PUBLIC "CGAL_LINKED_WITH_EMBREE")
#target_include_directories(example_advanced_features PUBLIC "${EMBREE_INCLUDE_DIRS}")
#target_link_libraries(example_advanced_features PUBLIC ${EMBREE_LIBRARIES})
  target_link_libraries(example_advanced_features PUBLIC CGAL::Embree_support)
else()
  message("NOTICE: Embree was not found. Embree ray intersection won't be available.")
endif()

#get_cmake_property(_variableNames VARIABLES)
#list (SORT _variableNames)
#foreach (_variableName ${_variableNames})
#    message(STATUS "${_variableName}=${${_variableName}}")
#endforeach()


find_package(PCL CONFIG)
if(PCL_FOUND)
    include_directories(${PCL_INCLUDE_DIRS})
    add_definitions(${PCL_DEFINITIONS})
    link_directories(${PCL_LIBRARY_DIRS})
    target_link_libraries(example_multi_dim_feature PUBLIC ${PCL_LIBRARIES})
else()
  message("NOTICE: PCL was not found. The PCL feature descriptor example won't be available.")
endif()

