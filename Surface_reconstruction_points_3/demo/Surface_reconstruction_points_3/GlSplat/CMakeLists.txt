
project(GlSplat)

cmake_minimum_required(VERSION 2.6.0)
if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}.${CMAKE_PATCH_VERSION}" VERSION_GREATER 2.8.3)
  cmake_policy(VERSION 2.8.4)
else()
  cmake_policy(VERSION 2.6)
endif()
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

set(QT_USE_QTOPENGL TRUE)
set(QT_USE_QTXML TRUE)

#New for Qt5 version !
if (NOT QT4 AND NOT QT5)
option (USE_QT5 "Use Qt5 version instead of Qt4" ON)
endif(NOT QT4 AND NOT QT5)

if(USE_QT5)
	
	find_package(Qt5 REQUIRED)

else (USE_QT5)

	find_package(Qt4 REQUIRED)

	if(CGAL_Qt4_FOUND AND QT4_FOUND)
		include(${QT_USE_FILE})
		message("Qt4 found")
		set(QT4 TRUE)
	endif()

endif (USE_QT5)

include(QtChoice)

find_package(GLEW REQUIRED)

include_directories( ${QT_LIBRARIES} ${GLEW_INCLUDE_DIR})

set(srcs Shader.cpp GlSplat.cpp)

qt_automoc( ${srcs})
qt_add_resources(srcs glsplat.qrc)

add_library(GlSplat SHARED ${srcs})
target_link_libraries(GlSplat ${QT_LIBRARIES} ${GLEW_LIBRARIES})

find_package(QGLViewer)
if(QGLVIEWER_FOUND)
  include_directories( ${QGLViewer_INCLUDES})
  add_executable(demo demo.cpp)
  target_link_libraries(demo ${QT_LIBRARIES} ${QGLViewer_LIBRARIES} GlSplat)
else()
  message("QGLViewer not found, demo won't be built")
endif()
