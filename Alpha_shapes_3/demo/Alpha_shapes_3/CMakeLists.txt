# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

project (Alpha_shape_3)
# Find includes in corresponding build directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)
# Instruct CMake to run moc automatically when needed.
set(CMAKE_AUTOMOC ON)
cmake_minimum_required(VERSION 2.6.2)
if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}.${CMAKE_PATCH_VERSION}" VERSION_GREATER 2.8.3)
  cmake_policy(VERSION 2.8.4)
else()
  cmake_policy(VERSION 2.6)
endif()

find_package(CGAL REQUIRED COMPONENTS Qt5)
include(${CGAL_USE_FILE})

find_package(Qt5 COMPONENTS Xml Script OpenGL)

find_package(OpenGL)
find_package(QGLViewer)

if ( CGAL_FOUND AND CGAL_Qt5_FOUND AND Qt5_FOUND AND OPENGL_FOUND AND QGLVIEWER_FOUND )


  include_directories (${QGLVIEWER_INCLUDE_DIR})
  include_directories (BEFORE ../../include ./ )

  # ui file, created wih Qt Designer
  qt5_wrap_ui( uis MainWindow.ui )

  # qrc files (resources files, that contain icons, at least)
  qt5_add_resources ( RESOURCE_FILES ./Alpha_shape_3.qrc )

  #qt4_automoc( MainWindow.cpp Viewer.cpp)

  add_executable  ( Alpha_shape_3 Alpha_shape_3.cpp MainWindow.cpp Viewer.cpp ${uis}  ${RESOURCE_FILES} )
  qt5_use_modules(Alpha_shape_3 Xml Script OpenGL)
  add_to_cached_list( CGAL_EXECUTABLE_TARGETS Alpha_shape_3 )

  target_link_libraries( Alpha_shape_3 ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES})
  target_link_libraries( Alpha_shape_3 ${QT_LIBRARIES} ${QGLVIEWER_LIBRARIES} )
  target_link_libraries( Alpha_shape_3 ${OPENGL_gl_LIBRARY} ${OPENGL_glu_LIBRARY} )

else()

  message(STATUS "NOTICE: This demo requires CGAL, the QGLViewer, OpenGL and Qt5, and will not be compiled.")

endif()
