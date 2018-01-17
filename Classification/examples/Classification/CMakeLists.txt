# Created by the script cgal_create_CMakeLists
# This is the CMake script for compiling a set of CGAL applications.

project( Classification_Examples )


cmake_minimum_required(VERSION 2.8.11)

# CGAL and its components
find_package( CGAL QUIET COMPONENTS  )

if ( NOT CGAL_FOUND )
  message(STATUS "This project requires the CGAL library, and will not be compiled.")
  return()  
endif()

# include helper file
include( ${CGAL_USE_FILE} )


# Boost and its components
find_package( Boost REQUIRED )

if ( NOT Boost_FOUND )
  message(STATUS "This project requires the Boost library, and will not be compiled.")
  return()  
endif()

find_package( Boost OPTIONAL_COMPONENTS serialization iostreams )

if( WIN32 )
# to avoid a warning with old cmake 
  set(_Boost_BZIP2_HEADERS             "boost/iostreams/filter/bzip2.hpp")
  set(_Boost_ZLIB_HEADERS              "boost/iostreams/filter/zlib.hpp")
  find_package( Boost OPTIONAL_COMPONENTS zlib)
endif()

find_package(TBB QUIET)

find_package(OpenCV QUIET)

# include for local directory
include_directories( BEFORE include )

# include for local package
include_directories( BEFORE ../../include )


# Creating entries for all C++ files with "main" routine
# ##########################################################

include( CGAL_CreateSingleSourceCGALProgram )

# Classification examples
set(targets
  example_classification
  example_ethz_random_forest
  example_feature
  example_generation_and_training)

# Classification requires some C++11 features
set(needed_cxx_features cxx_rvalue_references cxx_variadic_templates)

# Libraries and flags
set(classification_linked_libraries)
set(classification_compile_definitions)


if (Boost_SERIALIZATION_FOUND AND Boost_IOSTREAMS_FOUND)
  set(classification_linked_libraries ${classification_linked_libraries}
    ${Boost_SERIALIZATION_LIBRARY}
    ${Boost_IOSTREAMS_LIBRARY})
else()
  message(STATUS "NOTICE: Boost Serialization or IO Streams not found, no example will be compiled.")
  return()
endif()

if( WIN32 )
  if (Boost_ZLIB_FOUND)
    set(classification_linked_libraries ${classification_linked_libraries}
      ${Boost_ZLIB_LIBRARY})
  else()
    message(STATUS "NOTICE: Boost ZLIB not found, no example will be compiled.")
    return()
  endif()
endif()

find_package(OpenCV QUIET)
if (OpenCV_FOUND)
  message(STATUS "Found OpenCV ${OpenCV_VERSION}")
  include_directories(${OpenCV_INCLUDE_DIRS})
  set(classification_linked_libraries ${classification_linked_libraries}
    ${OpenCV_LIBS})
  set(classification_compile_definitions ${classification_compile_definitions}
    "-DCGAL_LINKED_WITH_OPENCV")

  set(targets ${targets} example_opencv_random_forest)
else()
  message(STATUS "NOTICE: OpenCV was not found. OpenCV random forest predicate for classification won't be available.")
endif()

# Creating targets with correct libraries and flags
foreach(target ${targets})
  create_single_source_cgal_program( "${target}.cpp" CXX_FEATURES ${needed_cxx_features} )
  if(TARGET ${target})
    target_link_libraries(${target} PUBLIC ${classification_linked_libraries})
    target_compile_definitions(${target} PUBLIC ${classification_compile_definitions})
    if(TBB_FOUND)
      CGAL_target_use_TBB(${target})
    endif()
  endif()
endforeach()





