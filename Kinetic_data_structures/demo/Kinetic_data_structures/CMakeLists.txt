# This is the CMake script for compiling a CGAL application.

project( Kinetic_data_structures_demo )

cmake_minimum_required(VERSION 2.8.11)
if(POLICY CMP0043)
  cmake_policy(SET CMP0043 OLD)
endif()

find_package(CGAL COMPONENTS Qt3 Core)

if ( CGAL_FOUND )

  include( ${CGAL_USE_FILE} )

  include( CGAL_CreateSingleSourceCGALProgram )
  
  find_package(Qt3-patched )
  # FindQt3-patched.cmake is FindQt3.cmake patched by CGAL developers, so
  # that it can be used together with Qt5: all its variables are prefixed
  # by "QT3_" instead of "QT_".
  
  if(CGAL_Qt3_FOUND AND QT3_FOUND)
  
    include( Qt3Macros-patched )
    qt3_automoc(  KDS_Delaunay_triangulation_2.cpp Delaunay_triangulation_3.cpp KDS_Delaunay_triangulation_stable_subset_2.cpp KDS_generate_data.cpp KDS_gui_2.cpp gui_3.cpp regular_triangulation_3.cpp SoQt_examiner_viewer.cpp )

    # Make sure the compiler can find generated .moc files
    include_directories(BEFORE ${CMAKE_CURRENT_BINARY_DIR})
  
    include_directories( ${QT3_INCLUDE_DIR} )

    include_directories (BEFORE include)
    include_directories (BEFORE ../../include)
    
    # Link with Boost.ProgramOptions (optional)
    find_package(Boost COMPONENTS program_options)
    if(Boost_PROGRAM_OPTIONS_FOUND)
      add_definitions( "-DCGAL_USE_BOOST_PROGRAM_OPTIONS" )
      list(APPEND CGAL_3RD_PARTY_LIBRARIES ${Boost_PROGRAM_OPTIONS_LIBRARY})
    endif()

    create_single_source_cgal_program( "KDS_Delaunay_triangulation_2.cpp" )
    create_single_source_cgal_program( "KDS_gui_2.cpp" )
    create_single_source_cgal_program( "KDS_generate_data.cpp")
    create_single_source_cgal_program( "KDS_Delaunay_triangulation_stable_subset_2.cpp" )
 
  else()
  
    message(STATUS "NOTICE: This demo requires Qt3 and the CGAL Qt3 library, and will not be compiled.")
  
  endif()
  
else()
  
    message(STATUS "NOTICE: This demo requires the CGAL library, and will not be compiled.")
  
endif()

