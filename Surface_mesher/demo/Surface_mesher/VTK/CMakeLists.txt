 cmake_minimum_required(VERSION 2.8.11)

PROJECT( mesh_a_3D_image )

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

set(PACKAGE_ROOT ../../..)

# Add several CGAL packages to the include and link paths,
# if they lie in ${PACKAGE_ROOT}/.
foreach(INC_DIR ${PACKAGE_ROOT}/include ${PACKAGE_ROOT}/../Mesh_2/include ${PACKAGE_ROOT}/../Data_structure_for_queries_3/include ${PACKAGE_ROOT}/../Marching_cube/include ${PACKAGE_ROOT}/../CGAL_ImageIO/include)
  if (EXISTS ${INC_DIR})
    include_directories (BEFORE ${INC_DIR})
  endif()
endforeach()
foreach(LIB_DIR ${PACKAGE_ROOT}/../CGAL_ImageIO/src/CGAL_ImageIO)
  if (EXISTS ${LIB_DIR})
    link_directories (${LIB_DIR})
  endif()
endforeach()

set( QT_USE_QTOPENGL TRUE )


find_package(VTK 6.0 COMPONENTS
  vtkRenderingOpenGL
  vtkGUISupportQt
  vtkIOImage
  vtkInteractionImage
  )
if(VTK_FOUND)
  add_definitions(-DCGAL_USE_VTK)
  include(${VTK_USE_FILE})
endif(VTK_FOUND)

if("${VTK_QT_VERSION}" STREQUAL "")
  message(FATAL_ERROR "VTK was not built with Qt")
endif()

if(VTK_QT_VERSION VERSION_GREATER "4")
  find_package(Qt5 COMPONENTS Core REQUIRED QUIET)

  find_package(CGAL REQUIRED ImageIO Qt5)
  IF(CGAL_FOUND AND CGAL_ImageIO_FOUND AND CGAL_Qt5_FOUND)

    add_definitions(-DQT_NO_KEYWORDS)
    include( ${CGAL_USE_FILE} )

    find_package(Qt5 QUIET COMPONENTS OpenGL Xml )
    find_package(OpenGL)

    SET (SRCS  mesh_a_3D_image.cpp)
    ADD_EXECUTABLE( mesh_a_3D_image MACOSX_BUNDLE ${SRCS})
    qt5_use_modules(mesh_a_3D_image Core Gui Widgets)
    TARGET_LINK_LIBRARIES( mesh_a_3D_image
      ${VTK_LIBRARIES}
      ${QT_LIBRARIES}
      vtkRenderingOpenGL
      vtkGUISupportQt
      vtkInteractionImage
      )

    SET (SRCS_VTK  mesh_a_VTK_3D_image.cpp)
    ADD_EXECUTABLE( mesh_a_VTK_3D_image MACOSX_BUNDLE ${SRCS_VTK})
    qt5_use_modules(mesh_a_VTK_3D_image Core Gui Widgets)
    TARGET_LINK_LIBRARIES( mesh_a_VTK_3D_image
      ${VTK_LIBRARIES}
      ${QT_LIBRARIES}
      vtkRenderingOpenGL
      vtkGUISupportQt
      vtkIOImage
      vtkInteractionImage
      )

  ENDIF(CGAL_FOUND AND CGAL_ImageIO_FOUND AND CGAL_Qt5_FOUND)

endif() #VTK_QT_VERSION VERSION_GREATER

