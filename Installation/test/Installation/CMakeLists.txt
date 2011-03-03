# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.


project( Installation_test )

CMAKE_MINIMUM_REQUIRED(VERSION 2.4.5)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)
 
if ( COMMAND cmake_policy )
  cmake_policy( SET CMP0003 NEW )  
endif()

macro(create_link_to_program COMPONENT )
  
  add_executable(link_to_${COMPONENT} link_to_${COMPONENT}.cpp )
  
  include_directories ( ${${COMPONENT}_3RD_PARTY_INCLUDE_DIRS} )     
        
  add_definitions( ${${COMPONENT}_3RD_PARTY_DEFINITIONS} )

  link_directories(  ${${COMPONENT}_3RD_PARTY_LIBRARIES_DIRS})
  
  # Link the executable to CGAL and third-party libraries
  if ( CGAL_AUTO_LINK_ENABLED )    
    target_link_libraries(link_to_${COMPONENT} ${CGAL_3RD_PARTY_LIBRARIES} ${${COMPONENT}_3RD_PARTY_LIBRARIES} )
  else()
    target_link_libraries(link_to_${COMPONENT} ${${COMPONENT}_LIBRARY} ${CGAL_3RD_PARTY_LIBRARIES} ${${COMPONENT}_3RD_PARTY_LIBRARIES} )
  endif()
  
endmacro()

find_package(CGAL QUIET)

if ( CGAL_FOUND )

  include( ${CGAL_USE_FILE} )

  include( CGAL_CreateSingleSourceCGALProgram )

  create_single_source_cgal_program( "endian.cpp" )

  create_link_to_program(CGAL)
 
  if ( WITH_CGAL_Core )
    create_link_to_program(CGAL_Core)
  endif()
  
  if ( WITH_CGAL_ImageIO )
    create_link_to_program(CGAL_ImageIO)
  endif()
  
  if ( WITH_CGAL_Qt3 )
    create_link_to_program(CGAL_Qt3)
  endif()
  
  if ( WITH_CGAL_Qt4 )
    find_package(Qt4 QUIET)
    if( QT4_FOUND )
      include(${QT_USE_FILE})
      create_link_to_program(CGAL_Qt4)
    endif()  
  endif()

  if(WIN32) 
    add_executable(test_gmp_mpfr_dll test_gmp_mpfr_dll.cpp)
    target_link_libraries(test_gmp_mpfr_dll Version)
  endif(WIN32)
  
else()
  
  message(STATUS "NOTICE: This program requires the CGAL library, and will not be compiled.")
  
endif()
