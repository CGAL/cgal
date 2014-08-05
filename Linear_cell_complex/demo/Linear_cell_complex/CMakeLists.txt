# Created by the script cgal_create_cmake_script (and then adapted manually).
# This is the CMake script for compiling a CGAL application.
# cmake ../ -DCMAKE_BUILD_TYPE=Debug

project (Linear_cell_complex_3_demo)

cmake_minimum_required(VERSION 2.6.2)
if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}.${CMAKE_PATCH_VERSION}" VERSION_GREATER 2.8.3)
  cmake_policy(VERSION 2.8.4)
else()
  cmake_policy(VERSION 2.6)
endif()

## To add expensive tests
# add_definitions("-DCGAL_CHECK_EXPENSIVE")

## For profilling with gprof
# add_definitions("-pg")
# SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pg")

## To add an is_valid test after each operation (only in debug mode)
# add_definitions("-DCGAL_CMAP_TEST_VALID_REMOVALS")
# add_definitions("-DCGAL_CMAP_TEST_VALID_CONTRACTIONS")
# add_definitions("-DCGAL_CMAP_TEST_VALID_INSERTIONS")

# Option allowing to profile each operation of the demo (cout times on stdout).
add_definitions(-DCGAL_PROFILE_LCC_DEMO)

##################

set( QT_USE_QTXML    TRUE )
set( QT_USE_QTMAIN   TRUE )
set( QT_USE_QTSCRIPT  TRUE )
set( QT_USE_QTOPENGL  TRUE )

find_package(CGAL COMPONENTS Qt4or5 REQUIRED)
include(${CGAL_USE_FILE})

find_package(Qt4or5)

find_package(OpenGL)
find_package(QGLViewer)

if ( NOT CGAL_FOUND )

  MESSAGE(STATUS "NOTICE: This demo requires CGAL "
    "and will not be compiled.")
elseif ( NOT QT${QT_VERSION_USED}_FOUND )
	
  MESSAGE(STATUS "NOTICE: This demo requires "
    "${QT_VERSION}, and will not be compiled.")
elseif ( NOT OPENGL_FOUND )

  MESSAGE(STATUS "NOTICE: This demo requires OpenGL "
    "and will not be compiled.")

elseif ( NOT QGLVIEWER_FOUND )

  MESSAGE(STATUS "NOTICE: This demo requires QGLViewer "
    " and will not be compiled. TODO correct this error")

else()
include_directories(${QGLVIEWER_INCLUDE_DIR})
include_directories(BEFORE . ../../include/)

FILE (GLOB CPP *.cpp)

# ui file, created wih Qt Designer

qt_wrap_ui(uis MainWindow.ui CreateMesh.ui CreateMenger.ui)

qt_add_resources (RESOURCE_FILES ./Linear_cell_complex_3.qrc)

qt_automoc(MainWindow.cpp Viewer.cpp)

add_executable(Linear_cell_complex_3_demo
			   ${CPP}
			   ${uis} 
			   ${RESOURCE_FILES} )

add_to_cached_list(CGAL_EXECUTABLE_TARGETS Linear_cell_complex_3_demo)

target_link_libraries(Linear_cell_complex_3_demo ${CGAL_LIBRARIES}
                                                 ${CGAL_3RD_PARTY_LIBRARIES})
target_link_libraries(Linear_cell_complex_3_demo ${QT_LIBRARIES}
                                                 ${QGLVIEWER_LIBRARIES} )
target_link_libraries(Linear_cell_complex_3_demo ${OPENGL_gl_LIBRARY}
                                                 ${OPENGL_glu_LIBRARY} )
endif()
