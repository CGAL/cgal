# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

project (Periodic_Lloyd_3)

# Find includes in corresponding build directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)
# Instruct CMake to run moc automatically when needed.
set(CMAKE_AUTOMOC ON)
cmake_minimum_required(VERSION 2.8.11)
if(POLICY CMP0043)
  cmake_policy(SET CMP0043 OLD)
endif()

find_package(CGAL COMPONENTS Qt5)
include(${CGAL_USE_FILE})


find_package(Qt5 QUIET COMPONENTS Xml Script Help OpenGL Svg)

find_package(OpenGL)

if(Qt5_FOUND)
  find_package(QGLViewer)
endif(Qt5_FOUND)

if ( CGAL_FOUND AND CGAL_Qt5_FOUND AND Qt5_FOUND AND OPENGL_FOUND AND QGLVIEWER_FOUND AND TARGET Qt5::qcollectiongenerator )


  include_directories (${QGLVIEWER_INCLUDE_DIR})
  include_directories (BEFORE ../../include ./ )

  # ui file, created wih Qt Designer
  qt5_wrap_ui( uis MainWindow.ui )

  # qrc files (resources files, that contain icons, at least)
  qt5_add_resources ( RESOURCE_FILES ./Periodic_Lloyd_3.qrc )


  if(DEFINED QT_QCOLLECTIONGENERATOR_EXECUTABLE)
  else()
    set(QT_QCOLLECTIONGENERATOR_EXECUTABLE qcollectiongenerator)
  endif()

  # generate QtAssistant collection file
  add_custom_command ( OUTPUT Periodic_Lloyd_3.qhc
      DEPENDS Periodic_Lloyd_3.qhp Periodic_Lloyd_3.qhcp
      COMMAND Qt5::qcollectiongenerator
                ${CMAKE_CURRENT_SOURCE_DIR}/Periodic_Lloyd_3.qhcp
                -o ${CMAKE_CURRENT_BINARY_DIR}/Periodic_Lloyd_3.qhc
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  )

  add_executable  ( Periodic_Lloyd_3 Periodic_Lloyd_3.qhc Periodic_Lloyd_3.cpp MainWindow.cpp Viewer.cpp ${uis}  ${RESOURCE_FILES} )

  qt5_use_modules(Periodic_Lloyd_3 Xml Script Help OpenGL Svg)
  add_to_cached_list( CGAL_EXECUTABLE_TARGETS Periodic_Lloyd_3 )

  target_link_libraries( Periodic_Lloyd_3 ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES})
  target_link_libraries( Periodic_Lloyd_3 ${QT_LIBRARIES} ${QGLVIEWER_LIBRARIES} )
  target_link_libraries( Periodic_Lloyd_3 ${OPENGL_gl_LIBRARY} ${OPENGL_glu_LIBRARY} )

else( CGAL_FOUND AND CGAL_Qt5_FOUND AND Qt5_FOUND AND OPENGL_FOUND AND QGLVIEWER_FOUND AND QT_QCOLLECTIONGENERATOR_EXECUTABLE )

  set(PERIODIC_LLOYD_MISSING_DEPS "")

  if(NOT CGAL_Qt5_FOUND)
    set(PERIODIC_LLOYD_MISSING_DEPS "the CGAL Qt5 library, ${PERIODIC_LLOYD_MISSING_DEPS}")
  endif()      	   

  if(NOT Qt5_FOUND)
    set(PERIODIC_LLOYD_MISSING_DEPS "Qt5, ${PERIODIC_LLOYD_MISSING_DEPS}")
  endif()      	   

  if(NOT OPENGL_FOUND)
    set(PERIODIC_LLOYD_MISSING_DEPS "OpenGL, ${PERIODIC_LLOYD_MISSING_DEPS}")
  endif()      	   

  if(NOT QGLVIEWER_FOUND)
    set(PERIODIC_LLOYD_MISSING_DEPS "QGLViewer, ${PERIODIC_LLOYD_MISSING_DEPS}")
  endif()      	   

  if(NOT QT_QCOLLECTIONGENERATOR_EXECUTABLE)
    set(PERIODIC_LLOYD_MISSING_DEPS "qcollectiongenerator, ${PERIODIC_LLOYD_MISSING_DEPS}")
  endif()


  message(STATUS "NOTICE: This demo requires ${PERIODIC_LLOYD_MISSING_DEPS}and will not be compiled.")

endif()
