cmake_minimum_required(VERSION 3.1...3.15)
project(Cone_spanners_2_Examples)

find_package(CGAL REQUIRED QUIET OPTIONAL_COMPONENTS Core)
include(${CGAL_USE_FILE})
find_package(LEDA QUIET)

if(CGAL_Core_FOUND OR LEDA_FOUND)
  if(MSVC)
    # Turn off a boost related warning that appears with VC2015
    # boost_1_65_1\boost\graph\named_function_params.hpp(240) :
    # warning C4172: returning address of local variable or temporary
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /wd4172")
  endif()

  # create a target per cppfile
  file(
    GLOB cppfiles
    RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)
  foreach(cppfile ${cppfiles})
    create_single_source_cgal_program("${cppfile}")
  endforeach()

else()

  message(
    STATUS
      "This program requires the CGAL_Core library (or LEDA), and will not be compiled."
  )

endif()
