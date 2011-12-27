project(image_to_vtk_viewer)

cmake_minimum_required(VERSION 2.6.2)
if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}" VERSION_GREATER 2.6)
  if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}.${CMAKE_PATCH_VERSION}" VERSION_GREATER 2.8.3)
    cmake_policy(VERSION 2.8.4)
  else()
    cmake_policy(VERSION 2.6)
  endif()
endif()


set(PACKAGE_ROOT ../..)

# Add several CGAL packages to the include and link paths,
# if they lie in ${PACKAGE_ROOT}/.
foreach(INC_DIR ${PACKAGE_ROOT}/include 
		${PACKAGE_ROOT}/../CGALimageIO/include )
  if(EXISTS ${INC_DIR})
    include_directories (BEFORE ${INC_DIR})
  endif()
endforeach()
foreach(LIB_DIR ${PACKAGE_ROOT}/../CGALimageIO/src/CGALImageIO)
  if (EXISTS ${LIB_DIR})
    link_directories (${LIB_DIR})
  endif()
endforeach()


find_package(CGAL REQUIRED ImageIO)
include( ${CGAL_USE_FILE} )
find_package(VTK QUIET)
find_package(Qt4 QUIET)

if(QT_FOUND AND VTK_FOUND)
  add_definitions(-DCGAL_USE_VTK)
  include(${VTK_USE_FILE})
  include(${QT_USE_FILE})

  add_definitions(${QT_DEFINITIONS})

  if(VTK_USE_QVTK)
    include_directories( ${VTK_QT_INCLUDE_DIR} )
    include_directories( ${QT_INCLUDE_DIR} )
    add_executable( image_to_vtk_viewer image_to_vtk_viewer.cpp )
    add_to_cached_list( CGAL_EXECUTABLE_TARGETS image_to_vtk_viewer )

    target_link_libraries( image_to_vtk_viewer
      QVTK
      vtkRendering
      vtkGraphics
      vtkIO
      vtkCommon
      ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES}
      ${VTK_QT_QT_LIBRARY}
      ${QT_LIBRARIES}
      )
  else(VTK_USE_QVTK)
    message(STATUS "NOTICE: This demo needs QVTK, and will not be compiled.")
  endif(VTK_USE_QVTK)
else()
  if(NOT VTK_FOUND)
    message(STATUS "NOTICE: This demo needs VTK, and will not be compiled.")
  endif()
  if(NOT QT_FOUND)
    message(STATUS "NOTICE: This demo needs Qt4, and will not be compiled.")
  endif()
endif()
