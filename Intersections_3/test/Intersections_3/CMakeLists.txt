# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.


project( Intersections_3_Tests )

cmake_minimum_required(VERSION 2.8.10)

find_package(CGAL QUIET)

if ( CGAL_FOUND )

  include( ${CGAL_USE_FILE} )

  include( CGAL_CreateSingleSourceCGALProgram )

  include_directories (BEFORE "../../include")

  # create a target per cppfile
  file(GLOB cppfiles RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)
  foreach(cppfile ${cppfiles})
    create_single_source_cgal_program( "${cppfile}" )
  endforeach()

foreach(VAL RANGE 1 55)
     message(STATUS "NOTICE: Compiling one_${VAL}")
    add_executable(one_${VAL} one.cpp)
    target_compile_definitions(one_${VAL} PRIVATE -DSEL=${VAL})

    target_link_libraries(one_${VAL} PRIVATE CGAL::CGAL)
    if(CGAL_3RD_PARTY_LIBRARIES)
      target_link_libraries(one_${VAL} PRIVATE ${CGAL_3RD_PARTY_LIBRARIES})
    endif()    
  endforeach()
  
else()
  
    message(STATUS "This program requires the CGAL library, and will not be compiled.")
  
endif()

