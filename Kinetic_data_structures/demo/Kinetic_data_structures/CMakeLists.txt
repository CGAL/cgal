# This is the CMake script for compiling a CGAL application.

project( Kinetic_data_structures_demo )

CMAKE_MINIMUM_REQUIRED(VERSION 2.4.5)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)
 
find_package(CGAL QUIET COMPONENTS Qt3 )

if ( CGAL_FOUND )

  include( ${CGAL_USE_FILE} )
  
  find_package(Qt3-patched QUIET )
  # FindQt3-patched.cmake is FindQt3.cmake patched by CGAL developers, so
  # that it can be used together with FindQt4: all its variables are prefixed
  # by "QT3_" instead of "QT_".
  
  if(CGAL_Qt3_FOUND AND QT3_FOUND)
  
    include( Qt3Macros-patched )
    qt3_automoc(  Delaunay_triangulation_2.cpp Delaunay_triangulation_3.cpp Delaunay_triangulation_stable_subset_2.cpp generate_data.cpp gui_2.cpp gui_3.cpp regular_triangulation_3.cpp SoQt_examiner_viewer.cpp )

    # Make sure the compiler can find generated .moc files
    include_directories(BEFORE ${CMAKE_CURRENT_BINARY_DIR})
  
    include_directories( ${QT3_INCLUDE_DIR} )

    include_directories (BEFORE include)

    find_package(Boost QUIET COMPONENTS program_options)

    if(Boost_PROGRAM_OPTIONS_FOUND)
      add_definitions(-DCGAL_USE_BOOST_PROGRAM_OPTIONS)
      include_directories(${BOOST_INCLUDE_DIRS})
    else(Boost_PROGRAM_OPTIONS_FOUND)
      message(STATUS "Warning: this directory uses Boost.Program_options. Without it, several features are disabled.")
    endif(Boost_PROGRAM_OPTIONS_FOUND)

    # Demo Delaunay_triangulation_2
    add_executable  (Delaunay_triangulation_2 Delaunay_triangulation_2.cpp)
    target_link_libraries(Delaunay_triangulation_2 ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES})
    if(Boost_PROGRAM_OPTIONS_FOUND)
      target_link_libraries(Delaunay_triangulation_2 ${Boost_PROGRAM_OPTIONS_LIBRARY})
    endif(Boost_PROGRAM_OPTIONS_FOUND)

    add_to_cached_list( CGAL_EXECUTABLE_TARGETS Delaunay_triangulation_2 )

    # Demo gui_2
    add_executable  (gui_2 gui_2.cpp)
    target_link_libraries(gui_2 ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES})
  
    add_to_cached_list( CGAL_EXECUTABLE_TARGETS gui_2 )
  
    # Demo generate_data
    add_executable  (generate_data generate_data.cpp)
    target_link_libraries(generate_data ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES})
    if(Boost_PROGRAM_OPTIONS_FOUND)
      target_link_libraries(generate_data ${Boost_PROGRAM_OPTIONS_LIBRARY})
    endif(Boost_PROGRAM_OPTIONS_FOUND)
  
    # Demo Delaunay_triangulation_stable_subset_2
    add_executable  (Delaunay_triangulation_stable_subset_2 Delaunay_triangulation_stable_subset_2.cpp)
    target_link_libraries(Delaunay_triangulation_stable_subset_2 ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES})
    if(Boost_PROGRAM_OPTIONS_FOUND)
      target_link_libraries(Delaunay_triangulation_stable_subset_2 ${Boost_PROGRAM_OPTIONS_LIBRARY})
    endif(Boost_PROGRAM_OPTIONS_FOUND)

    add_to_cached_list( CGAL_EXECUTABLE_TARGETS Delaunay_triangulation_stable_subset_2 )
  
  else()
  
    message(STATUS "NOTICE: This demo requires Qt3 and the CGAL Qt3 library, and will not be compiled.")
  
  endif()
  
else()
  
    message(STATUS "NOTICE: This demo requires the CGAL library, and will not be compiled.")
  
endif()

