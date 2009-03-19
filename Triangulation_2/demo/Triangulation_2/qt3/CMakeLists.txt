# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.


project( qt3_demo )

CMAKE_MINIMUM_REQUIRED(VERSION 2.4.5)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)
 
if ( COMMAND cmake_policy )
  cmake_policy( SET CMP0003 NEW )  
endif()
 
find_package(CGAL QUIET COMPONENTS Core Qt3 )

if ( CGAL_FOUND )

  include( ${CGAL_USE_FILE} )
  
  find_package(Qt3-patched QUIET )
  # FindQt3-patched.cmake is FindQt3.cmake patched by CGAL developers, so
  # that it can be used together with FindQt4: all its variables are prefixed
  # by "QT3_" instead of "QT_".
  
  if(CGAL_Qt3_FOUND AND QT3_FOUND)
  
    include( Qt3Macros-patched )
    qt3_automoc(  constrained.cpp constrained_delaunay_triangulation_2.cpp constrained_delaunay_triangulation_2_toolbar.cpp constrained_delaunay_triangulation_2_toolbar_layers.cpp delaunay_triangulation_2.cpp delaunay_triangulation_2_toolbar.cpp delaunay_triangulation_2_toolbar_layers.cpp regular_triangulation_2.cpp regular_triangulation_2_toolbar.cpp regular_triangulation_2_toolbar_layers.cpp triangulation_2_edit_vertex.cpp )

    # Make sure the compiler can find generated .moc files
    include_directories(BEFORE ${CMAKE_CURRENT_BINARY_DIR})
  
    include_directories( ${QT3_INCLUDE_DIR} )

    # Link the executables to CGAL and third-party libraries
    link_libraries(${QT3_LIBRARIES} ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} )

    add_executable(delaunay_triangulation_2 
      delaunay_triangulation_2.cpp
      triangulation_2_edit_vertex.cpp 
      delaunay_triangulation_2_toolbar.cpp 
      delaunay_triangulation_2_toolbar_layers.cpp)
  
    add_executable(regular_triangulation_2
      regular_triangulation_2.cpp
      triangulation_2_edit_vertex.cpp
      regular_triangulation_2_toolbar.cpp
      regular_triangulation_2_toolbar_layers.cpp)
    
    add_executable(constrained_delaunay_triangulation_2
      constrained_delaunay_triangulation_2.cpp 
      constrained_delaunay_triangulation_2_toolbar.cpp 
      constrained_delaunay_triangulation_2_toolbar_layers.cpp)
              
    add_executable(constrained constrained.cpp)
  
  else()
  
    message(STATUS "NOTICE: This demo requires Qt3 and the CGAL Qt3 library, and will not be compiled.")
  
  endif()
  
else()
  
    message(STATUS "NOTICE: This demo requires the CGAL library, and will not be compiled.")
  
endif()

