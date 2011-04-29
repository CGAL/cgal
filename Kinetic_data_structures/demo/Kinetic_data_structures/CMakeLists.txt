# This is the CMake script for compiling a CGAL application.

project( Kinetic_data_structures_demo )

cmake_minimum_required(VERSION 2.6.2)
if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}" VERSION_GREATER 2.6)
  if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}.${CMAKE_PATCH_VERSION}" VERSION_GREATER 2.8.3)
    cmake_policy(VERSION 2.8.4)
  else()
    cmake_policy(VERSION 2.6)
  endif()
endif()

find_package(CGAL COMPONENTS Qt3 Core)

if ( CGAL_FOUND )

  include( ${CGAL_USE_FILE} )
  
  find_package(Qt3-patched )
  # FindQt3-patched.cmake is FindQt3.cmake patched by CGAL developers, so
  # that it can be used together with FindQt4: all its variables are prefixed
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
      link_libraries( ${Boost_PROGRAM_OPTIONS_LIBRARY} )
    endif()

    # Demo KDS_Delaunay_triangulation_2
    add_executable  (KDS_Delaunay_triangulation_2 KDS_Delaunay_triangulation_2.cpp)
    target_link_libraries(KDS_Delaunay_triangulation_2 ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES})

    add_to_cached_list( CGAL_EXECUTABLE_TARGETS KDS_Delaunay_triangulation_2 )

    # Demo KDS_gui_2
    add_executable  (KDS_gui_2 KDS_gui_2.cpp)
    target_link_libraries(KDS_gui_2 ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES})
  
    add_to_cached_list( CGAL_EXECUTABLE_TARGETS KDS_gui_2 )
  
    # Demo KDS_generate_data
    add_executable  (KDS_generate_data KDS_generate_data.cpp)
    target_link_libraries(KDS_generate_data ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES})
    add_to_cached_list( CGAL_EXECUTABLE_TARGETS KDS_generate_data )

    # Demo KDS_Delaunay_triangulation_stable_subset_2
    add_executable  (KDS_Delaunay_triangulation_stable_subset_2 KDS_Delaunay_triangulation_stable_subset_2.cpp)
    target_link_libraries(KDS_Delaunay_triangulation_stable_subset_2 ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES})

    add_to_cached_list( CGAL_EXECUTABLE_TARGETS KDS_Delaunay_triangulation_stable_subset_2 )
  
  else()
  
    message(STATUS "NOTICE: This demo requires Qt3 and the CGAL Qt3 library, and will not be compiled.")
  
  endif()
  
else()
  
    message(STATUS "NOTICE: This demo requires the CGAL library, and will not be compiled.")
  
endif()

