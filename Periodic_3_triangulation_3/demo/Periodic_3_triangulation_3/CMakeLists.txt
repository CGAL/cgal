# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

project( Periodic_3_triangulation_3_demo ) 

CMAKE_MINIMUM_REQUIRED(VERSION 2.4.5)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

if(COMMAND cmake_policy)
  cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)

# Find CGAL
find_package(CGAL)
include( ${CGAL_USE_FILE} )

# Find Qt4 itself
set ( QT_USE_QTXML    TRUE )
set ( QT_USE_QTMAIN   TRUE )
set ( QT_USE_QTOPENGL TRUE )
set ( QT_USE_QTASSISTANT TRUE )
find_package(Qt4)

# Find OpenGL
find_package(OpenGL)

# Find QGLViewer
if(QT4_FOUND)
  include( ${QT_USE_FILE} )
  find_package(QGLViewer)
endif(QT4_FOUND)

if (CGAL_FOUND AND OPENGL_FOUND AND QT4_FOUND AND QGLVIEWER_FOUND) 

   include_directories( ${QGLVIEWER_INCLUDE_DIR} )

  # UI files (Qt Designer files)
  qt4_wrap_ui ( UI_FILES MainWindow.ui )

  # qrc files (resource files)
  qt4_add_resources( RESOURCE_FILES ./Periodic_3_triangulation_3.qrc )

  # use the Qt MOC preprocessor on classes that derive from QObject
  qt4_generate_moc ( Scene.h moc_Scene.cpp )
  qt4_generate_moc ( MainWindow.h moc_MainWindow.cpp )

  # Make sure the compiler can find generated .moc files
  include_directories(BEFORE ${CMAKE_CURRENT_BINARY_DIR})
  include_directories( ${QT_INCLUDE_DIR} )
  include_directories( ../../include )

  # The executable itself
  add_executable  ( periodic_3_triangulation_3_demo
  		    Scene.cpp moc_Scene.cpp
  		    periodic_3_triangulation_3_demo.cpp
  		    MainWindow.ui moc_MainWindow.cpp 
		    ${UI_FILES} ${RESOURCE_FILES})

  add_to_cached_list( CGAL_EXECUTABLE_TARGETS periodic_3_triangulation_3_demo)

  # Link the executable to CGAL and third-party libraries
  target_link_libraries(periodic_3_triangulation_3_demo ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} ${QT_LIBRARIES} ${QGLVIEWER_LIBRARIES} ${OPENGL_gl_LIBRARY} ${OPENGL_glu_LIBRARY} )

else (CGAL_FOUND AND OPENGL_FOUND AND QT4_FOUND AND QGLVIEWER_FOUND) 

  set(PERIODIC_TRIANGULATION_MISSING_DEPS "")

  if(NOT CGAL_FOUND)
    set(PERIODIC_TRIANGULATION_MISSING_DEPS "the CGAL Qt4 library, ${PERIODIC_TRIANGULATION_MISSING_DEPS}")
  endif()

  if(NOT QT4_FOUND)
    set(PERIODIC_TRIANGULATION_MISSING_DEPS "Qt4, ${PERIODIC_TRIANGULATION_MISSING_DEPS}")
  endif()

  if(NOT OPENGL_FOUND)
    set(PERIODIC_TRIANGULATION_MISSING_DEPS "OpenGL, ${PERIODIC_TRIANGULATION_MISSING_DEPS}")
  endif()

  if(NOT QGLVIEWER_FOUND)
    set(PERIODIC_TRIANGULATION_MISSING_DEPS "QGLViewer, ${PERIODIC_TRIANGULATION_MISSING_DEPS}")
  endif()

  message(STATUS "NOTICE: This demo requires ${PERIODIC_TRIANGULATION_MISSING_DEPS}and will not be compiled.")

endif (CGAL_FOUND AND OPENGL_FOUND AND QT4_FOUND AND QGLVIEWER_FOUND) 
