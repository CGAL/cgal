project(image_to_vtk_viewer)

cmake_minimum_required(VERSION 2.8.11)

find_package(CGAL REQUIRED ImageIO)
include( ${CGAL_USE_FILE} )
find_package(VTK QUIET NO_MODULE)

if(VTK_FOUND)
  add_definitions(-DCGAL_USE_VTK)

  if(TARGET vtkRenderingQt AND TARGET vtkFiltersModeling)
    find_package(VTK COMPONENTS vtkRenderingQt vtkFiltersModeling NO_MODULE)
    include(${VTK_USE_FILE})
    find_package(Qt${VTK_QT_VERSION} COMPONENTS QtGui)
    if(NOT TARGET Qt${VTK_QT_VERSION}::QtGui)
      message(STATUS "NOTICE: vtkRenderingQt needs Qt${VTK_QT_VERSION}, and will not be compiled.")
      return()
    endif()

    add_executable( image_to_vtk_viewer image_to_vtk_viewer.cpp )
    add_to_cached_list( CGAL_EXECUTABLE_TARGETS image_to_vtk_viewer )

    target_link_libraries( image_to_vtk_viewer
      ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES}
      ${VTK_LIBRARIES} Qt${VTK_QT_VERSION}::QtGui
      )
  else()
    message(STATUS "NOTICE: This demo needs vtkRenderingQt and vtkFiltersModeling, and will not be compiled.")
  endif()
else()
  message(STATUS "NOTICE: This demo needs VTK, and will not be compiled.")
endif()
