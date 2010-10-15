PROJECT( mesh_a_3D_image )

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

set(PACKAGE_ROOT ../../..)

# Add several CGAL packages to the include and link paths,
# if they lie in ${PACKAGE_ROOT}/.
foreach(INC_DIR ${PACKAGE_ROOT}/include ${PACKAGE_ROOT}/../Mesh_2/include ${PACKAGE_ROOT}/../Data_structure_for_queries_3/include ${PACKAGE_ROOT}/../Marching_cube/include ${PACKAGE_ROOT}/../CGALimageIO/include)
  if (EXISTS ${INC_DIR})
    include_directories (BEFORE ${INC_DIR})
  endif()
endforeach()
foreach(LIB_DIR ${PACKAGE_ROOT}/../CGALimageIO/src/CGALimageIO)
  if (EXISTS ${LIB_DIR})
    link_directories (${LIB_DIR})
  endif()
endforeach()


FIND_PACKAGE(CGAL)
IF(CGAL_FOUND)
  include( ${CGAL_USE_FILE} )
ENDIF(CGAL_FOUND)

FIND_PACKAGE(VTK)
ADD_DEFINITIONS(-DCGAL_USE_VTK)
IF(NOT VTK_DIR)
  MESSAGE(FATAL_ERROR "Please set VTK_DIR.")
ENDIF(NOT VTK_DIR)
INCLUDE(${VTK_USE_FILE})

SET(QT_QMAKE_EXECUTABLE ${VTK_QT_QMAKE_EXECUTABLE} CACHE FILEPATH "")
SET(QT_MOC_EXECUTABLE ${VTK_QT_MOC_EXECUTABLE} CACHE FILEPATH "")
SET(QT_UIC_EXECUTABLE ${VTK_QT_UIC_EXECUTABLE} CACHE FILEPATH "")
SET(DESIRED_QT_VERSION ${VTK_DESIRED_QT_VERSION} CACHE FILEPATH "")
FIND_PACKAGE(Qt)
IF(QT_USE_FILE)
  INCLUDE(${QT_USE_FILE})
ELSE(QT_USE_FILE)
  SET(QT_LIBRARIES   ${QT_QT_LIBRARY})
ENDIF(QT_USE_FILE)

SET (SRCS
  mesh_a_3D_image.cpp
)

SET (SRCS_VTK
  mesh_a_VTK_3D_image.cpp
)

# Use the include path and library for Qt that is used by VTK.
INCLUDE_DIRECTORIES( ${QT_INCLUDE_DIR} ${QT_QTGUI_INCLUDE_DIR}
                     ${QT_QTCORE_INCLUDE_DIR})

ADD_EXECUTABLE( mesh_a_3D_image MACOSX_BUNDLE ${SRCS})

TARGET_LINK_LIBRARIES( mesh_a_3D_image
  QVTK
  ${QT_LIBRARIES}
  vtkRendering
  vtkGraphics
  vtkIO
  vtkCommon
  CGAL
  CGAL-ImageIO
)

SET (SRCS_VTK
  mesh_a_VTK_3D_image.cpp
)

ADD_EXECUTABLE( mesh_a_VTK_3D_image MACOSX_BUNDLE ${SRCS_VTK})

TARGET_LINK_LIBRARIES( mesh_a_VTK_3D_image
  QVTK
  ${QT_LIBRARIES}
  vtkRendering
  vtkGraphics
  vtkIO
  vtkCommon
  CGAL
  CGAL-ImageIO
)



