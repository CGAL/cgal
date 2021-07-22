cmake_minimum_required(VERSION 3.1...3.20)
project(CGALimageIO_Demo)

if(NOT POLICY CMP0070 AND POLICY CMP0053)
  # Only set CMP0053 to OLD with CMake<3.10, otherwise there is a warning.
  cmake_policy(SET CMP0053 OLD)
endif()

set(PACKAGE_ROOT ../..)

# Add several CGAL packages to the include and link paths,
# if they lie in ${PACKAGE_ROOT}/.
foreach(INC_DIR ${PACKAGE_ROOT}/include ${PACKAGE_ROOT}/../CGAL_ImageIO/include)
  if(EXISTS ${INC_DIR})
    include_directories(BEFORE ${INC_DIR})
  endif()
endforeach()
foreach(LIB_DIR ${PACKAGE_ROOT}/../CGAL_ImageIO/src/CGAL_ImageIO)
  if(EXISTS ${LIB_DIR})
    link_directories(${LIB_DIR})
  endif()
endforeach()

find_package(CGAL REQUIRED ImageIO)
find_package(VTK QUIET NO_MODULE)

if(VTK_FOUND)
  add_definitions(-DCGAL_USE_VTK)

  if(TARGET vtkRenderingQt AND TARGET vtkFiltersModeling)
    find_package(VTK COMPONENTS vtkRenderingQt vtkFiltersModeling NO_MODULE)
    include(${VTK_USE_FILE})
    find_package(Qt${VTK_QT_VERSION} COMPONENTS Gui)
    if(NOT TARGET Qt${VTK_QT_VERSION}::Gui)
      message(
        STATUS
          "NOTICE: vtkRenderingQt needs Qt${VTK_QT_VERSION}, and will not be compiled."
      )
      return()
    endif()

    add_executable(image_to_vtk_viewer image_to_vtk_viewer.cpp)
    add_to_cached_list(CGAL_EXECUTABLE_TARGETS image_to_vtk_viewer)

    target_link_libraries(
      image_to_vtk_viewer ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES}
      ${VTK_LIBRARIES} Qt${VTK_QT_VERSION}::Gui)
  else()
    message(
      STATUS
        "NOTICE: This demo needs vtkRenderingQt and vtkFiltersModeling, and will not be compiled."
    )
  endif()
else()
  message(STATUS "NOTICE: This demo needs VTK, and will not be compiled.")
endif()
