# CMake
# Fernando de Goes (fdegoes@caltech.edu)
# Copyright @ 2011
project(GSoC_pwsrec2D)

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

#########
# Files #
#########

#set( 
#HDRS 
#console_color.h 
#Cost.h 
#Dynamic_priority_queue.h 
#Reconstruction_triangulation_2.h
#Reconstruction_edge_2.h 
#Reconstruction_simplification_2.h 
#Reconstruction_simplification_kerneled_2.h 
#Reconstruction_vertex_base_2.h
#Reconstruction_face_base_2.h
#Sample.h

#)

set(
SRCS 
test_basic.cpp
)

#########
# Build #
#########

# Includes
include_directories(../../include/CGAL)
include_directories(BEFORE . ./build)
include(${QT_USE_FILE})
include_directories(${GLUT_INCLUDE_DIR})
include_directories(${OPENGL_INCLUDE_DIR})
include(${CGAL_USE_FILE})


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
