# Created by the script cgal_create_cmake_script (and then adapted manually).
# This is the CMake script for compiling a CGAL application.


project( mesh_2_qt3_demo )

cmake_minimum_required(VERSION 2.6.2)
if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}" VERSION_GREATER 2.6)
  if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}.${CMAKE_PATCH_VERSION}" VERSION_GREATER 2.8.3)
    cmake_policy(VERSION 2.8.4)
  else()
    cmake_policy(VERSION 2.6)
  endif()
endif()

find_package(CGAL COMPONENTS Core Qt3 )
include( ${CGAL_USE_FILE} )

find_package(Qt3-patched )
# FindQt3-patched.cmake is FindQt3.cmake patched by CGAL developers, so
# that it can be used together with FindQt4: all its variables are prefixed
# by "QT3_" instead of "QT_".

if ( CGAL_FOUND AND QT3_FOUND AND CGAL_Qt3_FOUND)

  include( Qt3Macros-patched )
  QT3_AUTOMOC(  icons.cpp mesh_2_demo.cpp Qt_widget_styled_layer.cpp Qt_widget_style_editor.cpp Show_points.cpp Show_segments_base.cpp )

  # Make sure the compiler can find generated .moc files
  include_directories(BEFORE ${CMAKE_CURRENT_BINARY_DIR})

  include_directories( ${QT3_INCLUDE_DIR} )

  add_executable  (mesh_2_demo icons.cpp mesh_2_demo.cpp Qt_widget_styled_layer.cpp Qt_widget_style_editor.cpp Show_points.cpp Show_segments_base.cpp)

  add_to_cached_list( CGAL_EXECUTABLE_TARGETS mesh_2_demo )

  # Link the main executable to CGAL and third-party libraries
  target_link_libraries(mesh_2_demo ${QT3_LIBRARIES} ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} )

else()

  message(STATUS "NOTICE: This demo requires CGAL and Qt3, and will not be compiled.")

endif()
