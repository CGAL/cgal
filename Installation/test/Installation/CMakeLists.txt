# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.


project( Installation_Tests )

cmake_minimum_required(VERSION 2.8.11)


macro(create_link_to_program COMPONENT )
  
  add_executable(link_to_${COMPONENT} link_to_${COMPONENT}.cpp )
  
  include_directories ( ${${COMPONENT}_3RD_PARTY_INCLUDE_DIRS} )     
        
  add_definitions( ${${COMPONENT}_3RD_PARTY_DEFINITIONS} )

  link_directories(  ${${COMPONENT}_3RD_PARTY_LIBRARIES_DIRS})
  
  # Link the executable to CGAL and third-party libraries
  if ( CGAL_AUTO_LINK_ENABLED OR CGAL_HEADER_ONLY)
    target_link_libraries(link_to_${COMPONENT} ${CGAL_3RD_PARTY_LIBRARIES} ${${COMPONENT}_3RD_PARTY_LIBRARIES} )
  else()
    target_link_libraries(link_to_${COMPONENT} CGAL::${COMPONENT} ${CGAL_3RD_PARTY_LIBRARIES} ${${COMPONENT}_3RD_PARTY_LIBRARIES} )
  endif()
  
endmacro()

find_package(CGAL QUIET COMPONENTS Core)

if ( CGAL_FOUND )

  include( ${CGAL_USE_FILE} )

  include( CGAL_CreateSingleSourceCGALProgram )

  create_single_source_cgal_program( "endian.cpp" )

  create_single_source_cgal_program( "deprecation_warning.cpp" )

  create_single_source_cgal_program( "test_use_h.cpp" )

  create_single_source_cgal_program( will_fail.cpp will_fail_aux.cpp )
  if(BUILD_TESTING)
    set_property(TEST execution___of__will_fail PROPERTY WILL_FAIL TRUE)
  endif()

  find_package( TBB QUIET )
  create_single_source_cgal_program( "test_TBB.cpp" )
  if(TBB_FOUND)
    CGAL_target_use_TBB(test_TBB)
  endif()

  create_link_to_program(CGAL)
 
  if ( WITH_CGAL_Core )
    create_link_to_program(CGAL_Core)
  endif()
  
  if ( WITH_CGAL_ImageIO )
    create_link_to_program(CGAL_ImageIO)
  endif()

  if ( WITH_CGAL_Qt4 )
    find_package(Qt4 QUIET)
    if( QT4_FOUND )
      include(${QT_USE_FILE})
      create_link_to_program(CGAL_Qt4)
    endif()  
  endif()

  if ( WITH_CGAL_Qt5 )
    find_package(Qt5 QUIET)
    if( Qt5_FOUND )
      create_link_to_program(CGAL_Qt5)
      if (CGAL_HEADER_ONLY)
        qt5_use_modules(link_to_CGAL_Qt5)
      endif()
    endif()
  endif()

  if(WIN32) 
    add_executable(test_gmp_mpfr_dll test_gmp_mpfr_dll.cpp)
    target_link_libraries(test_gmp_mpfr_dll Version)
  endif(WIN32)

  find_package( LEDA QUIET)
  if(LEDA_FOUND)
    include( ${LEDA_USE_FILE} )
    add_executable(test_LEDA_IO test_LEDA_IO.cpp)
    target_link_libraries(test_LEDA_IO ${LEDA_LIBRARIES})
  endif()

else()
  
  message(STATUS "NOTICE: This program requires the CGAL library, and will not be compiled.")
  
endif()
