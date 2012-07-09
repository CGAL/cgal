# Created by the script cgal_create_cmake_script (and then adapted manually).
# This is the CMake script for compiling a CGAL application.
# cmake ../ -DCMAKE_BUILD_TYPE=Debug

project (Linear_cell_complex_3_demo)

cmake_minimum_required(VERSION 2.6.2)
if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}" VERSION_GREATER 2.6)
  if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}.${CMAKE_PATCH_VERSION}" VERSION_GREATER 2.8.3)
    cmake_policy(VERSION 2.8.4)
  else()
    cmake_policy(VERSION 2.6)
  endif()
endif()

find_package(CGAL COMPONENTS Qt4)
include(${CGAL_USE_FILE})

set( QT_USE_QTXML    TRUE )
set( QT_USE_QTMAIN   TRUE )
set( QT_USE_QTSCRIPT  TRUE )
set( QT_USE_QTOPENGL  TRUE )
find_package(Qt4)

find_package(OpenGL)
find_package(QGLViewer)

if ( NOT (CGAL_FOUND AND CGAL_Qt4_FOUND AND QT4_FOUND AND 
      OPENGL_FOUND AND QGLVIEWER_FOUND) )

  MESSAGE(STATUS "NOTICE: This demo requires CGAL, QGLViewer, OpenGL "
    "and Qt4, and will not be compiled.")

else()

include(${QT_USE_FILE})
include_directories(${QGLVIEWER_INCLUDE_DIR})
include_directories(BEFORE . ../../include/)

# Option allowing to profile each operation of the demo.
add_definitions(-DCGAL_PROFILE_LCC_DEMO)

# ui file, created wih Qt Designer
qt4_wrap_ui(uis MainWindow.ui CreateMesh.ui CreateMenger.ui)

# qrc files (resources files, that contain icons, at least)
qt4_add_resources (RESOURCE_FILES ./Linear_cell_complex_3.qrc)

qt4_automoc(MainWindow.cpp Viewer.cpp)

add_executable(Linear_cell_complex_3_demo
  Linear_cell_complex_3_demo.cpp MainWindow.cpp
  Viewer.cpp Linear_cell_complex_3_subdivision.cpp
  Linear_cell_complex_pqq_subdivision.cpp
  ${uis} ${RESOURCE_FILES} )

add_to_cached_list(CGAL_EXECUTABLE_TARGETS Linear_cell_complex_3_demo)

target_link_libraries(Linear_cell_complex_3_demo ${CGAL_LIBRARIES}
                                                 ${CGAL_3RD_PARTY_LIBRARIES})
target_link_libraries(Linear_cell_complex_3_demo ${QT_LIBRARIES}
                                                 ${QGLVIEWER_LIBRARIES} )
target_link_libraries(Linear_cell_complex_3_demo ${OPENGL_gl_LIBRARY}
                                                 ${OPENGL_glu_LIBRARY} )

endif()
