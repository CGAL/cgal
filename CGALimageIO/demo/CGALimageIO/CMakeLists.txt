project(image_to_vtk_viewer)

cmake_minimum_required(VERSION 2.4.5)

if ( COMMAND cmake_policy )
  cmake_policy( SET CMP0003 NEW )  
endif()

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

  if(VTK_USE_QVTK)
    add_executable( image_to_vtk_viewer image_to_vtk_viewer.cpp)

    target_link_libraries( image_to_vtk_viewer
      QVTK
      vtkRendering
      vtkGraphics
      vtkIO
      vtkCommon
      ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES}
      )
  else(VTK_USE_QVTK)
    message(STATUS "NOTICE: This demo needs QVTK, and will not be compiled.")
  endif(VTK_USE_QVTK)
else(VTK_FOUND AND WITH_VTK)

  message(STATUS "NOTICE: This demo needs VTK, and will not be compiled.")

endif(VTK_FOUND AND WITH_VTK)
