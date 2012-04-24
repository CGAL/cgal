# This file must be included in your CMakeLists.txt to use
# cgal_map_viewer_vtk.h.
# You need to link the libraries in your executable by using 
# TARGET_LINK_LIBRARIES( myexec ${MAP_VIEWER_LIBRARIES})

if ( NOT CGAL_FOUND )
  MESSAGE(FATAL_ERROR "Please install CGAL.")
ENDIF(NOT CGAL_FOUND)

include( ${CGAL_USE_FILE} )

set(QT_USE_QT3SUPPORT true)
set(QT_USE_QTOPENGL true)
FIND_PACKAGE(Qt4 REQUIRED)

INCLUDE(${QT_USE_FILE})
ADD_DEFINITIONS(${QT_DEFINITIONS})

FIND_PACKAGE(VTK REQUIRED)
IF(NOT VTK_DIR)
  MESSAGE(FATAL_ERROR "Please set VTK_DIR.")
ENDIF(NOT VTK_DIR)

INCLUDE(${VTK_USE_FILE})

SET(QT_QMAKE_EXECUTABLE ${VTK_QT_QMAKE_EXECUTABLE} CACHE FILEPATH "")

# Use the include path and library for Qt that is used by VTK.
INCLUDE_DIRECTORIES(${QT_INCLUDE_DIR} ${CMAKE_CURRENT_BINARY_DIR}
                    ${CMAKE_CURRENT_SOURCE_DIR})

SET (MAP_VIEWER_LIBRARIES
     QVTK ${QT_LIBRARIES} vtkRendering vtkGraphics
     vtkIO vtkCommon)

ADD_DEFINITIONS("-DCGAL_LCC_USE_VIEWER  -DCGAL_LCC_USE_VTK")
    
message(STATUS "Libraries for map_viewer found. You need to link them "
        "in your executable by using "
        "TARGET_LINK_LIBRARIES( myexec \${MAP_VIEWER_LIBRARIES})")
