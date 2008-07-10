# This is the CMake script for compiling a CGAL application.

project( Polyhedron ) 

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

foreach(INCDIR ../../../GraphicsView/include ../../../Surface_mesher/include)
  if(EXISTS ${INCDIR})
    include_directories (BEFORE ${INCDIR})
  endif()
endforeach()

include_directories(include)

if(EXISTS ../../../GraphicsView/src/CGALQt4)
  link_directories(../../../GraphicsView/src/CGALQt4)
endif(EXISTS ../../../GraphicsView/src/CGALQt4)

find_package(CGAL REQUIRED)
include( ${CGAL_USE_FILE} )

set( QT_USE_QTXML    TRUE )
set( QT_USE_QTMAIN   TRUE )
set( QT_USE_QTSCRIPT  TRUE )
set( QT_USE_QTOPENGL  TRUE )

find_package(Qt4       REQUIRED)
include(${QT_USE_FILE})

find_package(QGLViewer REQUIRED )
include_directories( ${QGLVIEWER_INCLUDE_DIR} )
    
qt4_wrap_ui( UI_FILES MainWindow.ui )

include(AddFileDependencies)

qt4_generate_moc( MainWindow.h MainWindow_moc.cpp )
add_file_dependencies(  MainWindow_moc.cpp MainWindow.h )

qt4_generate_moc( Viewer.h Viewer_moc.cpp )
add_file_dependencies( Viewer_moc.cpp Viewer.h )

qt4_generate_moc( Scene.h Scene_moc.cpp )
add_file_dependencies( Scene_moc.cpp Scene.h )

qt4_add_resources ( RESOURCE_FILES Polyhedron_3.qrc )

add_executable  ( Polyhedron_3 MainWindow_moc.cpp Scene_moc.cpp Viewer_moc.cpp main.cpp Viewer.cpp ${UI_FILES} ${RESOURCE_FILES} )

# Link with Qt libraries
target_link_libraries( Polyhedron_3 ${QT_LIBRARIES} )

if ( NOT AUTO_LINK_ENABLED )    
  # If auto-link is not enabled (on Linux for example), link with libCGAL explicitly.
  target_link_libraries( Polyhedron_3 CGAL )
endif()
target_link_libraries( Polyhedron_3 CGAL-Qt4 ${QGLVIEWER_LIBRARY} ${OPENGL_gl_LIBRARY} ${OPENGL_glu_LIBRARY} )
