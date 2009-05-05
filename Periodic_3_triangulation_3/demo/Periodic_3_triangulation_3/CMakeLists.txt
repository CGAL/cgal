# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

project( Triangulation_3 ) 

CMAKE_MINIMUM_REQUIRED(VERSION 2.4.5)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

find_package(CGAL REQUIRED)
include( ${CGAL_USE_FILE} )

set ( QT_USE_QTXML    TRUE )
set ( QT_USE_QTMAIN   TRUE )
set ( QT_USE_QTOPENGL TRUE )

find_package(OpenGL REQUIRED)

find_package(Qt4 REQUIRED)
include( ${QT_USE_FILE} )

find_package(QGLViewer REQUIRED)
include_directories( ${QGLVIEWER_INCLUDE_DIR} )


qt4_wrap_ui ( UI_FILES periodic.ui )

qt4_generate_moc ( Scene.h moc_Scene.cpp )

# Make sure the compiler can find generated .moc files
include_directories(BEFORE ${CMAKE_CURRENT_BINARY_DIR})
include_directories( ${QT_INCLUDE_DIR} )
include_directories( ../../include )

add_executable  ( demo_delaunay_qt
		  Scene.cpp moc_Scene.cpp
		  demo_delaunay_qt.cpp
		  periodic.ui )

# TODO: what to do if AUTO_LINK_ENABLED?
# TODO: replace CGAL and CGAL-core++ by CGAL_[3RD_PARTY_]LIBRARIES
# Link the executable to CGAL and third-party libraries
if ( AUTO_LINK_ENABLED )
   #target_link_libraries(demo_delaunay_qt ${CGAL_3RD_PARTY_LIBRARIES} ${QT_LIBRARIES} )
else()
  target_link_libraries(demo_delaunay_qt ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} ${QT_LIBRARIES} ${QGLVIEWER_LIBRARIES} ${OPENGL_gl_LIBRARY} ${OPENGL_glu_LIBRARY} )
#  target_link_libraries(demo_delaunay_qt ${CGAL_LIBRARIES}
#  ${CGAL_3RD_PARTY_LIBRARIES} ${QT_LIBRARIES} ${QGLVIEWER_LIBRARIES}
#  ${OPENGL_gl_LIBRARY} ${OPENGL_glu_LIBRARY} )
endif()

