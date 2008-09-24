project(image_to_vtk_viewer)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

set(PACKAGE_ROOT ../..)

# Add several CGAL packages to the include and link paths,
# if they lie in ${PACKAGE_ROOT}/.
foreach(INC_DIR ${PACKAGE_ROOT}/include 
		${PACKAGE_ROOT}/../CGALimageIO/include )
  if(EXISTS ${INC_DIR})
    include_directories (BEFORE ${INC_DIR})
  endif()
endforeach()
foreach(LIB_DIR ${PACKAGE_ROOT}/../CGALimageIO/src/CGALimageIO)
  if (EXISTS ${LIB_DIR})
    link_directories (${LIB_DIR})
  endif()
endforeach()


find_package(CGAL REQUIRED ImageIO)
include( ${CGAL_USE_FILE} )
find_package(VTK QUIET)

if(VTK_FOUND AND WITH_VTK)
  add_definitions(-DCGAL_USE_VTK)
  include(${VTK_USE_FILE})

  set(QT_QMAKE_EXECUTABLE ${VTK_QT_QMAKE_EXECUTABLE} CACHE FILEPATH "")
  set(QT_MOC_EXECUTABLE ${VTK_QT_MOC_EXECUTABLE} CACHE FILEPATH "")
  set(QT_UIC_EXECUTABLE ${VTK_QT_UIC_EXECUTABLE} CACHE FILEPATH "")
  set(DESIRED_QT_VERSION ${VTK_DESIRED_QT_VERSION} CACHE FILEPATH "")

  find_package(Qt QUIET)
  if(QT_FOUND)
    include(${QT_USE_FILE})

    # Use the include path and library for Qt that is used by VTK.
    include_directories( ${QT_INCLUDE_DIR} ${QT_QTGUI_INCLUDE_DIR}
      ${QT_QTCORE_INCLUDE_DIR})

    add_executable( image_to_vtk_viewer image_to_vtk_viewer.cpp)

    target_link_libraries( image_to_vtk_viewer
      QVTK
      ${QT_LIBRARIES}
      vtkRendering
      vtkGraphics
      vtkIO
      vtkCommon
      ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES}
      )
  else(QT_FOUND)
    message(STATUS "NOTICE: This demo needs Qt, and will not be compiled.")
  endif(QT_FOUND)
else(VTK_FOUND)

  message(STATUS "NOTICE: This demo needs VTK, and will not be compiled.")

endif(VTK_FOUND)
