# CMake
# Fernando de Goes (fdegoes@caltech.edu)
# Copyright @ 2011
project(pwsrec2D)

cmake_minimum_required(VERSION 2.4.5)
cmake_policy(VERSION 2.4.5)

set( CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true )

############
# Packages #
############

set( QT_USE_QTXML    TRUE )
set( QT_USE_QTMAIN   TRUE )
set( QT_USE_QTSCRIPT  TRUE )
set( QT_USE_QTOPENGL  TRUE )

find_package(OpenGL)
find_package(Qt4)
find_package(CGAL COMPONENTS Qt4)



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

#########
# Build #
#########

# Includes
include_directories(BEFORE . ./build)
include(${QT_USE_FILE})
include_directories(${GLUT_INCLUDE_DIR})
include_directories(${OPENGL_INCLUDE_DIR})
include(${CGAL_USE_FILE})



  include_directories (BEFORE "../../include")
  include_directories (BEFORE "include")



# UI files
qt4_wrap_ui( DT_UI_FILES ${UIS} )

# QRC files
qt4_add_resources( DT_RESOURCE_FILES pwsrec.qrc )

# MOC files
qt4_generate_moc( window.h moc_window.cxx )
qt4_generate_moc( glviewer.h moc_glviewer.cxx )
qt4_generate_moc( dialog_options.h moc_dialog_options.cxx )

# The executable itself.
add_executable( ${PROJECT_NAME} ${SRCS} ${MOCS} ${DT_UI_FILES} ${DT_RESOURCE_FILES} )

# Link with Qt libraries
target_link_libraries( ${PROJECT_NAME} ${QT_LIBRARIES} )
    
# Link with Glut and OpenGL
target_link_libraries( ${PROJECT_NAME} ${GLUT_LIBRARY} ${OPENGL_LIBRARY} )

# Link with CImg dependencies
if( NOT WIN32 )
	target_link_libraries( ${PROJECT_NAME} -L/usr/X11R6/lib -lm -lpthread -lX11 )
endif()
    
# Link with CGAL
target_link_libraries( ${PROJECT_NAME} ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES})
