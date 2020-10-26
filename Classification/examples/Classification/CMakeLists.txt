# Created by the script cgal_create_CMakeLists
# This is the CMake script for compiling a set of CGAL applications.

cmake_minimum_required(VERSION 3.1...3.15)
project(Classification_Examples)

# CGAL and its components
find_package(CGAL REQUIRED)

# Boost and its components
find_package(Boost REQUIRED)

if(NOT Boost_FOUND)
  message(
    STATUS "This project requires the Boost library, and will not be compiled.")
  return()
endif()

set(Classification_dependencies_met TRUE)

find_package(Boost OPTIONAL_COMPONENTS serialization iostreams)
include(CGAL_Boost_serialization_support)
include(CGAL_Boost_iostreams_support)

if(NOT TARGET CGAL::Boost_serialization_support)
  message(
    STATUS
      "NOTICE: This project requires Boost Serialization, and will not be compiled."
  )
  set(Classification_dependencies_met FALSE)
endif()
if(NOT TARGET CGAL::Boost_iostreams_support)
  message(
    STATUS
      "NOTICE: This project requires Boost IO Streams, and will not be compiled."
  )
  set(Classification_dependencies_met FALSE)
endif()

find_package(OpenCV QUIET COMPONENTS core ml) # Need core + machine learning
include(CGAL_OpenCV_support)
if(NOT TARGET CGAL::OpenCV_support)
  message(
    STATUS
      "NOTICE: OpenCV was not found. OpenCV random forest predicate for classification won't be available."
  )
endif()

find_package(TensorFlow QUIET)
include(CGAL_TensorFlow_support)
if(NOT TARGET CGAL::TensorFlow_support)
  message(
    STATUS
      "NOTICE: TensorFlow was not found. TensorFlow neural network predicate for classification won't be available."
  )
endif()

find_package(Eigen3 3.1.0 REQUIRED) #(3.1.0 or greater)
include(CGAL_Eigen_support)
if(NOT TARGET CGAL::Eigen_support)
  message(
    STATUS "This project requires the Eigen library, and will not be compiled.")
  set(Classification_dependencies_met FALSE)
endif()

find_package(TBB QUIET)
include(CGAL_TBB_support)

if(NOT Classification_dependencies_met)
  return()
endif()

create_single_source_cgal_program( "example_classification.cpp" )
create_single_source_cgal_program( "example_ethz_random_forest.cpp" )
create_single_source_cgal_program( "example_feature.cpp" )
create_single_source_cgal_program( "example_generation_and_training.cpp" )
create_single_source_cgal_program( "example_mesh_classification.cpp" )
create_single_source_cgal_program( "example_cluster_classification.cpp" )
create_single_source_cgal_program( "gis_tutorial_example.cpp" )
create_single_source_cgal_program( "example_deprecated_conversion.cpp" )

if (TARGET CGAL::OpenCV_support)
  create_single_source_cgal_program( "example_opencv_random_forest.cpp" )
  target_link_libraries(example_opencv_random_forest
                        PUBLIC CGAL::OpenCV_support)
endif()

if(TARGET CGAL::TensorFlow_support)
  create_single_source_cgal_program("example_tensorflow_neural_network.cpp")
  target_link_libraries(example_opencv_random_forest
                        PUBLIC CGAL::TensorFlow_support)
endif()


foreach(target
    example_classification
    example_ethz_random_forest
    example_feature
    example_generation_and_training
    example_mesh_classification
    example_cluster_classification
    example_opencv_random_forest
    example_tensorflow_neural_network
    gis_tutorial_example
    example_deprecated_conversion)
  if(TARGET ${target})
    target_link_libraries(
      ${target} PUBLIC CGAL::Eigen_support CGAL::Boost_iostreams_support
                       CGAL::Boost_serialization_support)
    if(TARGET CGAL::TBB_support)
      target_link_libraries(${target} PUBLIC CGAL::TBB_support)
    endif()
  endif()
endforeach()
