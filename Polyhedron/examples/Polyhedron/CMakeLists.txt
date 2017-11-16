# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.


project( Polyhedron_Examples )

cmake_minimum_required(VERSION 2.8.10)

# If you want to visualize a polyhedron, you can use the following viewer
# based on qt. Just uncomment the following line, plus the line find_package(CGALBasicViewer) below
# find_package(CGAL COMPONENTS Qt5)

# If you don't want to visualize, use the following line (otherwise comment it)
find_package(CGAL QUIET)

if ( CGAL_FOUND )

  include( ${CGAL_USE_FILE} )
  # find_package(CGALBasicViewer)

  include( CGAL_CreateSingleSourceCGALProgram )

  include_directories (BEFORE "../../include")

  
  # create a target per cppfile
  file(GLOB cppfiles RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)
  foreach(cppfile ${cppfiles})
    create_single_source_cgal_program( "${cppfile}" )

    if (USE_BASIC_VIEWER)
      get_filename_component(exe_name "${cppfile}" NAME_WE)
      target_link_libraries("${exe_name}" ${BASIC_VIEWER_LIBRARIES})
      qt5_use_modules("${exe_name}" ${BASIC_VIEWER_MODULES})
    endif(USE_BASIC_VIEWER)
    
  endforeach()

else()
  
    message(STATUS "This program requires the CGAL library, and will not be compiled.")
  
endif()

