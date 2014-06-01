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

<<<<<<< HEAD
#New for Qt5 version !
=======
>>>>>>> 1f669ad498d664482b6b486e3e95a50a56076759
option (USE_QT5 "Use Qt5 version instead of Qt4" OFF)

if(USE_QT5)
	if(WIN32)
	message("Qt5 on Windows needs Windows SDK.")
	
	set(CMAKE_LIBRARY_PATH "C:\\Program Files (x86)\\Windows Kits\\8.1\\Lib\\winv6.3\\um\\x64")
	
	endif()
<<<<<<< HEAD
	find_package(CGAL COMPONENTS Qt5)
=======
	find_package(CGAL)
>>>>>>> 1f669ad498d664482b6b486e3e95a50a56076759
	include(${CGAL_USE_FILE})

	FIND_PACKAGE(Qt5Core REQUIRED)
	FIND_PACKAGE(Qt5Gui REQUIRED)
	FIND_PACKAGE(Qt5OpenGL REQUIRED)
	FIND_PACKAGE(Qt5Widgets REQUIRED)
	FIND_PACKAGE(Qt5Xml REQUIRED)
<<<<<<< HEAD
	FIND_PACKAGE(Qt5Quick REQUIRED)
=======
	# FIND_PACKAGE(Qt5Quick REQUIRED)
	# FIND_PACKAGE(Qt5Xml REQUIRED)
>>>>>>> 1f669ad498d664482b6b486e3e95a50a56076759
	
	if (Qt5Core_FOUND AND Qt5Gui_FOUND AND Qt5OpenGL_FOUND AND Qt5Widgets_FOUND AND Qt5Xml_FOUND)
		message("Qt5 found")
		set(QT5 TRUE)
		set(CMAKE_AUTOMOC ON)
<<<<<<< HEAD

=======
		# list(APPEND QT_LIBRARIES ${Qt5Widgets_LIBRARIES} ${Qt5Quick_LIBRARIES} ${Qt5Xml_LIBRARIES})
>>>>>>> 1f669ad498d664482b6b486e3e95a50a56076759
list(APPEND QT_INCLUDE_DIR ${Qt5Core_INCLUDE_DIRS} ${Qt5Gui_INCLUDE_DIRS} ${Qt5Widgets_INCLUDE_DIRS} ${Qt5OpenGL_INCLUDE_DIRS} ${Qt5Xml_INCLUDE_DIRS}) 

list(APPEND QT_LIBRARIES ${Qt5Core_LIBRARIES} ${Qt5Gui_LIBRARIES} ${Qt5Widgets_LIBRARIES} ${Qt5OpenGL_LIBRARIES} ${Qt5Xml_LIBRARIES})

list(APPEND QT_DEFINITIONS ${Qt5Core_DEFINITIONS} ${Qt5Gui_DEFINITIONS} ${Qt5Widgets_DEFINITIONS} ${Qt5OpenGL_DEFINITIONS} ${Qt5Xml_DEFINITIONS})
	endif ()
<<<<<<< HEAD
else (USE_QT5)
=======
else ()
>>>>>>> 1f669ad498d664482b6b486e3e95a50a56076759
	find_package(CGAL COMPONENTS Qt4)

	include(${CGAL_USE_FILE})

	find_package(Qt4)
	include(${QT_USE_FILE})

	if(CGAL_Qt4_FOUND AND QT4_FOUND)
		message("Qt4 found")
		set(QT4 TRUE)
	endif()
endif (USE_QT5)


set( QT_USE_QTXML    TRUE )
set( QT_USE_QTMAIN   TRUE )
set( QT_USE_QTSCRIPT  TRUE )
set( QT_USE_QTOPENGL  TRUE )

find_package(OpenGL)
find_package(QGLViewer)

if ( NOT CGAL_FOUND )

  MESSAGE(STATUS "NOTICE: This demo requires CGAL "
    "and will not be compiled.")
elseif ( NOT (QT4 OR QT5) )
	
  MESSAGE(STATUS "NOTICE: This demo requires "
    "Qt4 or Qt5, and will not be compiled.")
elseif ( NOT OPENGL_FOUND )

  MESSAGE(STATUS "NOTICE: This demo requires OpenGL "
    "and will not be compiled.")
<<<<<<< HEAD
elseif ( NOT QGLVIEWER_FOUND )

  MESSAGE(STATUS "NOTICE: This demo requires QGLViewer "
    " and will not be compiled. TODO correct this error")
else ()

=======
elseif (QGL_VIEWER_FOUND )

  MESSAGE(STATUS "NOTICE: This demo requires QGLViewer "
    " and will not be compiled. TODO correct this error")
else ()
>>>>>>> 1f669ad498d664482b6b486e3e95a50a56076759
include_directories(${QGLVIEWER_INCLUDE_DIR})
include_directories(BEFORE . ../../include/)

FILE (GLOB CPP *.cpp)

if(QT4)
	# ui file, created wih Qt Designer
	qt4_wrap_ui(uis MainWindow.ui CreateMesh.ui CreateMenger.ui)

	# qrc files (resources files, that contain icons, at least)
	qt4_add_resources (RESOURCE_FILES ./Linear_cell_complex_3.qrc)

	qt4_automoc(MainWindow.cpp Viewer.cpp)
endif()
if(QT5)

	qt5_wrap_ui(uis MainWindow.ui CreateMesh.ui CreateMenger.ui)

	qt5_add_resources (RESOURCE_FILES ./Linear_cell_complex_3.qrc)

endif()
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
<<<<<<< HEAD

endif()
=======
endif()
>>>>>>> 1f669ad498d664482b6b486e3e95a50a56076759
