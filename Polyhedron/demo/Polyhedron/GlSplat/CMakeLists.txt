
project(GlSplat)

cmake_minimum_required(VERSION 2.6.0)
if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}.${CMAKE_PATCH_VERSION}" VERSION_GREATER 2.8.3)
  cmake_policy(VERSION 2.8.4)
else()
  cmake_policy(VERSION 2.6)
endif()
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

find_package(Qt4 REQUIRED)
find_package(GLEW REQUIRED)

set(QT_USE_QTOPENGL TRUE)
set(QT_USE_QTXML TRUE)
include(${QT_USE_FILE})
include_directories( ${QT_QTOPENGL_INCLUDE_DIR} ${QT_QTXML_INCLUDE_DIR} ${GLEW_INCLUDE_DIR})

set(srcs Shader.cpp GlSplat.cpp)

qt4_automoc( ${srcs})
qt4_add_resources(srcs glsplat.qrc)

add_library(GlSplat SHARED ${srcs})
target_link_libraries(GlSplat ${QT_QTCORE_LIBRARY} ${QT_QTGUI_LIBRARY} ${QT_QTOPENGL_LIBRARY} ${GLEW_LIBRARIES})

find_package(QGLViewer)
if(QGLVIEWER_FOUND)
  include_directories( ${QGLViewer_INCLUDES})
  add_executable(demo demo.cpp)
  target_link_libraries(demo ${QT_QTXML_LIBRARY} ${QGLViewer_LIBRARIES} GlSplat)
else()
  message("QGLViewer not found, demo won't be built")
endif()
