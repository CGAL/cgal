# Created by the script cgal_create_CMakeLists
# This is the CMake script for compiling a set of CGAL applications.

cmake_minimum_required(VERSION 3.1...3.15)
project( Classification_Examples )


# CGAL and its components
find_package( CGAL QUIET COMPONENTS  )

if ( NOT CGAL_FOUND )
  message(STATUS "This project requires the CGAL library, and will not be compiled.")
  return()
endif()

# Boost and its components
find_package( Boost REQUIRED )

if ( NOT Boost_FOUND )
  message(STATUS "This project requires the Boost library, and will not be compiled.")
  return()
endif()

set(Classification_dependencies_met TRUE)

find_package( Boost OPTIONAL_COMPONENTS serialization iostreams )
if (NOT Boost_SERIALIZATION_FOUND)
  message(STATUS "NOTICE: This project requires Boost Serialization, and will not be compiled.")
  set(Classification_dependencies_met FALSE)
endif()
if (NOT Boost_IOSTREAMS_FOUND)
  message(STATUS "NOTICE: This project requires Boost IO Streams, and will not be compiled.")
  set(Classification_dependencies_met FALSE)
endif()

find_package(OpenCV QUIET COMPONENTS core ml) # Need core + machine learning
if (NOT OpenCV_FOUND)
  message(STATUS "NOTICE: OpenCV was not found. OpenCV random forest predicate for classification won't be available.")
endif()

find_package(TensorFlow QUIET)
if (NOT TensorFlow_FOUND)
  message(STATUS "NOTICE: TensorFlow was not found. TensorFlow neural network predicate for classification won't be available.")
endif()

find_package(Eigen3 3.1.0 REQUIRED) #(3.1.0 or greater)
if (NOT EIGEN3_FOUND)
  message(STATUS "This project requires the Eigen library, and will not be compiled.")
  set(Classification_dependencies_met FALSE)
endif()

find_package(TBB QUIET)

if (NOT Classification_dependencies_met)
  return()
endif()

create_single_source_cgal_program( "example_classification.cpp" )
create_single_source_cgal_program( "example_ethz_random_forest.cpp" )
create_single_source_cgal_program( "example_feature.cpp" )
create_single_source_cgal_program( "example_generation_and_training.cpp" )
create_single_source_cgal_program( "example_mesh_classification.cpp" )
create_single_source_cgal_program( "example_cluster_classification.cpp" )
create_single_source_cgal_program( "example_deprecated_conversion.cpp" )

if (OpenCV_FOUND)
  create_single_source_cgal_program( "example_opencv_random_forest.cpp" )
  CGAL_target_use_OpenCV(example_opencv_random_forest)
endif()

if (TensorFlow_FOUND)
  create_single_source_cgal_program( "example_tensorflow_neural_network.cpp" )
  CGAL_target_use_TensorFlow(example_tensorflow_neural_network)
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
    example_deprecated_conversion)
  if(TARGET ${target})
    CGAL_target_use_Eigen(${target})
    CGAL_target_use_Boost_IOStreams(${target})
    CGAL_target_use_Boost_Serialization(${target})
    if(TBB_FOUND)
      CGAL_target_use_TBB(${target})
    endif()
  endif()
endforeach()
