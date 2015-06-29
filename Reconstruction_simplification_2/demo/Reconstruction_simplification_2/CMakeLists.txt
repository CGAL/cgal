# This is the CMake script for compiling the Reconstruction_simplification_2 demo.

project(Reconstruction_simplification_2_demo)

cmake_minimum_required(VERSION 2.6.2)
if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}.${CMAKE_PATCH_VERSION}" VERSION_GREATER 2.8.3)
  cmake_policy(VERSION 2.8.4)
else()
  cmake_policy(VERSION 2.6)
endif()

set( CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true )

# Include this package's headers first
include_directories( BEFORE ./ ./include ../../include )

# Find CGAL and CGAL Qt4
find_package(CGAL COMPONENTS Qt4)
include( ${CGAL_USE_FILE} )

# Find Qt4 itself
set( QT_USE_QTXML    TRUE )
set( QT_USE_QTMAIN   TRUE )
set( QT_USE_QTSCRIPT  TRUE )
set( QT_USE_QTOPENGL  TRUE )
find_package(Qt4)

# Find OpenGL
find_package(OpenGL)

# Find QGLViewer
if(QT4_FOUND)
  include(${QT_USE_FILE})
  find_package(QGLViewer )
endif(QT4_FOUND)

if(CGAL_FOUND AND CGAL_Qt4_FOUND AND QT4_FOUND AND OPENGL_FOUND AND QGLVIEWER_FOUND)

  include_directories ( ${QGLVIEWER_INCLUDE_DIR} )

  set(
  SRCS 
  glviewer.cpp 
  main.cpp 
  window.cpp 
  render.cpp 
  )

  set( 
  MOCS 
  moc_dialog_options.cxx 
  moc_glviewer.cxx 
  moc_window.cxx 
  )

  set(
  UIS
  pwsrec.ui 
  options.ui
  )
  
  qt4_wrap_ui( UI_FILES ${UIS} )

  include(AddFileDependencies)

  qt4_generate_moc( "window.h" "${CMAKE_CURRENT_BINARY_DIR}/moc_window.cxx" )
  add_file_dependencies( moc_window.cxx "${CMAKE_CURRENT_SOURCE_DIR}/window.h" )

  qt4_generate_moc( "glviewer.h" "${CMAKE_CURRENT_BINARY_DIR}/moc_glviewer.cxx" )
  add_file_dependencies( moc_glviewer.cxx "${CMAKE_CURRENT_SOURCE_DIR}/glviewer.h" )

  qt4_generate_moc( "dialog_options.h" "${CMAKE_CURRENT_BINARY_DIR}/moc_dialog_options.cxx" )
  add_file_dependencies( moc_dialog_options.cxx "${CMAKE_CURRENT_SOURCE_DIR}/dialog_options.h" )

  qt4_add_resources ( RESOURCE_FILES pwsrec.qrc )

  add_executable  ( Reconstruction_simplification_2_demo ${SRCS} ${MOCS} ${UI_FILES} ${RESOURCE_FILES} ${CGAL_ADDITIONAL_FILES} ${CGAL_RESOURCE_FILES})

  # Link with Qt libraries
  target_link_libraries( Reconstruction_simplification_2_demo ${QT_LIBRARIES} )

  # Link with CGAL
  target_link_libraries( Reconstruction_simplification_2_demo ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} )

  # Link with libQGLViewer, OpenGL
  target_link_libraries( Reconstruction_simplification_2_demo ${QGLVIEWER_LIBRARIES} ${OPENGL_gl_LIBRARY} ${OPENGL_glu_LIBRARY} )

  # Link with CImg dependencies
  if( NOT WIN32 )
	  target_link_libraries( Reconstruction_simplification_2_demo -L/usr/X11R6/lib -lm -lpthread -lX11 )
  endif()
    
  add_to_cached_list( CGAL_EXECUTABLE_TARGETS Reconstruction_simplification_2_demo )

else (CGAL_FOUND AND CGAL_CGAL_Qt4_FOUND AND QT4_FOUND AND OPENGL_FOUND AND QGLVIEWER_FOUND)

  set(RS2_MISSING_DEPS "")

  if(NOT CGAL_FOUND)
    set(RS2_MISSING_DEPS "the CGAL library, ${RS2_MISSING_DEPS}")
  endif()

  if(NOT CGAL_Qt4_FOUND)
    set(RS2_MISSING_DEPS "the CGAL Qt4 library, ${RS2_MISSING_DEPS}")
  endif()

  if(NOT QT4_FOUND)
    set(RS2_MISSING_DEPS "Qt4, ${RS2_MISSING_DEPS}")
  endif()

  if(NOT OPENGL_FOUND)
    set(RS2_MISSING_DEPS "OpenGL, ${RS2_MISSING_DEPS}")
  endif()

  if(NOT QGLVIEWER_FOUND)
    set(RS2_MISSING_DEPS "QGLViewer, ${RS2_MISSING_DEPS}")
  endif()

  message(STATUS "NOTICE: This demo requires ${RS2_MISSING_DEPS}and will not be compiled.")

endif (CGAL_FOUND AND CGAL_Qt4_FOUND AND QT4_FOUND AND OPENGL_FOUND AND QGLVIEWER_FOUND)
