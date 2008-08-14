# This is the CMake script for compiling a CGAL application.

project( Polyhedron ) 

cmake_minimum_required(VERSION 2.4)

option(ENABLE_FORWARD_DECL "Enable the forward declaration of Polyhedron and Nef_polyhedron." OFF)

if(COMMAND cmake_policy)
  cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

foreach(INCDIR ./include ../../../GraphicsView/include ../../../Surface_mesher/include )
  if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${INCDIR})
    include_directories (BEFORE ${CMAKE_CURRENT_SOURCE_DIR}/${INCDIR})
  endif()
endforeach()

if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/../../../GraphicsView/src/CGALQt4)
  link_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../../GraphicsView/src/CGALQt4)
endif(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/../../../GraphicsView/src/CGALQt4)

find_package(CGAL REQUIRED)
include( ${CGAL_USE_FILE} )

find_package(Boost 1.34.1 COMPONENTS thread REQUIRED)

set( QT_USE_QTXML    TRUE )
set( QT_USE_QTMAIN   TRUE )
set( QT_USE_QTSCRIPT  TRUE )
set( QT_USE_QTOPENGL  TRUE )

find_package(OpenGL REQUIRED)

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

if(ENABLE_FORWARD_DECL)
  add_definitions(-DUSE_FORWARD_DECL)
  add_executable  ( Polyhedron_3 MainWindow.cpp
                                 main.cpp
                                 Viewer.cpp
                                 Scene.cpp
                                 Scene_nef_polyhedron_operations.cpp
                                 Scene_polyhedron_operations.cpp
                                 Scene_nef_and_polyhedron_operations.cpp
				 Scene_rendering.cpp
				 Scene_nef_rendering.cpp
				 MainWindow_boolean_operations.cpp
				 MainWindow_convex_hull.cpp
				 MainWindow_curvature_estimation.cpp
				 MainWindow_inside_out.cpp
				 MainWindow_kernel.cpp
				 MainWindow_pca.cpp
				 MainWindow_remeshing.cpp
				 MainWindow_self_intersection.cpp
				 MainWindow_simplify.cpp
				 MainWindow_subdivision_methods.cpp

                                 MainWindow_moc.cpp
                                 Scene_moc.cpp
                                 Viewer_moc.cpp
				 texture.cpp

                               ${UI_FILES} ${RESOURCE_FILES} )
else(ENABLE_FORWARD_DECL)
  add_file_dependencies( main.cpp MainWindow_moc.cpp Scene_moc.cpp Viewer_moc.cpp )
  add_executable  ( Polyhedron_3 main.cpp ${UI_FILES} ${RESOURCE_FILES} )
endif(ENABLE_FORWARD_DECL)


# Link with Qt libraries
target_link_libraries( Polyhedron_3 ${QT_LIBRARIES} )

# If auto-link is not enabled (on Linux for example), link with libCGAL and Boost.Thread explicitly.
if ( NOT AUTO_LINK_ENABLED )    
  target_link_libraries( Polyhedron_3 CGAL ${Boost_THREAD_LIBRARY_RELEASE} )
endif()

# Link with libCGAL-Qt4, libQGLViewer, OpenGL
target_link_libraries( Polyhedron_3 CGAL-Qt4 ${QGLVIEWER_LIBRARY} ${OPENGL_gl_LIBRARY} ${OPENGL_glu_LIBRARY} )

# Blas and Lapack, under Linux
if( NOT AUTO_LINK_ENABLED )
  add_definitions( -DCGAL_USE_F2C )
  target_link_libraries( Polyhedron_3 lapack blas )
endif()
