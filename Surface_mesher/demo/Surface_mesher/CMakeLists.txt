set(prj Surface_mesher)

cmake_minimum_required(VERSION 3.1...3.15)
project(Surface_mesher_Demo)

if(NOT POLICY CMP0070 AND POLICY CMP0053)
  # Only set CMP0053 to OLD with CMake<3.10, otherwise there is a warning.
  cmake_policy(SET CMP0053 OLD)
endif()

# Find includes in corresponding build directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)
# Instruct CMake to run moc automatically when needed.
set(CMAKE_AUTOMOC ON)

if(POLICY CMP0071)
  cmake_policy(SET CMP0071 NEW)
endif()
if(POLICY CMP0072)
  # About the use of OpenGL
  cmake_policy(SET CMP0072 NEW)
endif()

if(POLICY CMP0074)
  cmake_policy(SET CMP0074 NEW)
endif()

set(PACKAGE_ROOT ${CMAKE_CURRENT_SOURCE_DIR}/../..)

# Add several CGAL packages to the include and link paths,
# if they lie in ${PACKAGE_ROOT}/.

foreach(LIB_DIR ${PACKAGE_ROOT}/../CGAL_ImageIO/src/CGAL_ImageIO
                ${PACKAGE_ROOT}/../GraphicsView/src/CGALQt5)
  if(EXISTS ${LIB_DIR})
    link_directories(${LIB_DIR})
  endif()
endforeach()

include_directories(./)

# QGLViwer needs Qt5 configured with QtOpenGL and QtXml support

find_package(CGAL REQUIRED COMPONENTS ImageIO Qt5)
if(CGAL_Qt5_FOUND AND CGAL_ImageIO_FOUND)

  add_definitions(-DQT_NO_KEYWORDS)

  find_package(Qt5 QUIET COMPONENTS OpenGL Xml Svg)
  find_package(OpenGL)

  if(Qt5_FOUND
     AND OPENGL_FOUND
     AND OPENGL_GLU_FOUND)

    set(sources
        Raw_image_dialog.cpp
        colorlisteditor.cpp
        values_list.cpp
        mainwindow.cpp
        Surface_mesher.cpp
        viewer.cpp
        volume.cpp)

    if(EXISTS ${PACKAGE_ROOT}/../Marching_cube)
      option(SURFACE_MESH_DEMO_USE_MARCHING_CUBE
             "Embed a marching cube implementation in the Surface Mesh demo."
             ON)
      mark_as_advanced(SURFACE_MESH_DEMO_USE_MARCHING_CUBE)
    endif()

    if(EXISTS ${PACKAGE_ROOT}/include/CGAL/Polyhedral_surface_3.h)
      option(SURFACE_MESH_DEMO_WITH_POLYHEDRAL_SURFACE
             "Compile the support for polyhedral surfaces." OFF)
      mark_as_advanced(SURFACE_MESH_DEMO_WITH_POLYHEDRAL_SURFACE)
    endif()

    option(
      SURFACE_MESH_DEMO_VERBOSE
      "Set this option if you want the Surface Mesh demo to display messages on standard output."
      OFF)
    mark_as_advanced(SURFACE_MESH_DEMO_VERBOSE)

    if(SURFACE_MESH_DEMO_VERBOSE)
      add_definitions(-DCGAL_SURFACE_MESHER_VERBOSE)
    endif()

    if(SURFACE_MESH_DEMO_WITH_POLYHEDRAL_SURFACE)
      set(sources ${sources} polyhedral_surface.cpp)
    else(SURFACE_MESH_DEMO_WITH_POLYHEDRAL_SURFACE)
      add_definitions(-DCGAL_DO_NOT_USE_POLYHEDRAL_SURFACE)
    endif(SURFACE_MESH_DEMO_WITH_POLYHEDRAL_SURFACE)

    if(SURFACE_MESH_DEMO_USE_MARCHING_CUBE)
      set(sources ${sources} ${PACKAGE_ROOT}/../Marching_cube/src/mc/ply.c)
      add_definitions(-DCGAL_SURFACE_MESH_DEMO_USE_MARCHING_CUBE)
    endif()

    qt5_generate_moc("surface.h" "${CMAKE_CURRENT_BINARY_DIR}/moc_surface.cpp")

    set(sources ${sources} moc_surface.cpp)

    qt5_wrap_ui(uis ui/values_list.ui ui/mainwindow.ui ui/optionsdialog.ui
                ui/raw_image.ui)

    qt5_add_resources(CGAL_Qt5_RESOURCE_FILES values_list.qrc
                      surface_mesher.qrc)

    add_executable(${prj} ${sources} ${uis} ${CGAL_Qt5_RESOURCE_FILES}
                          ${CGAL_Qt5_MOC_FILES})

    add_to_cached_list(CGAL_EXECUTABLE_TARGETS ${prj})

    find_package(VTK QUIET COMPONENTS vtkImagingGeneral vtkIOImage NO_MODULE)
    if(VTK_FOUND)
      if(VTK_USE_FILE)
        include(${VTK_USE_FILE})
      endif()
      if("${VTK_VERSION_MAJOR}" GREATER "5" OR VTK_VERSION VERSION_GREATER 5)
        message(STATUS "VTK found")
        add_definitions(-DCGAL_USE_VTK)
        if(TARGET VTK::IOImage)
          set(VTK_LIBRARIES VTK::IOImage VTK::ImagingGeneral)
        endif()
      else()
        message(STATUS "Vtk must be at least Rel 6")
      endif()
    else()
      message(STATUS "For reading Dicom files install VTK first")
    endif()

    target_link_libraries(
      ${prj} PRIVATE CGAL::CGAL CGAL::CGAL_Qt5 CGAL::CGAL_ImageIO
                     ${OPENGL_LIBRARIES} ${VTK_LIBRARIES})

    include(${CGAL_MODULES_DIR}/CGAL_add_test.cmake)
    cgal_add_compilation_test(${prj})

    include(${CGAL_MODULES_DIR}/CGAL_add_test.cmake)
    cgal_add_compilation_test(${prj})
  endif(
    Qt5_FOUND
    AND OPENGL_FOUND
    AND OPENGL_GLU_FOUND)
else(
  CGAL_Qt5_FOUND
  AND CGAL_ImageIO_FOUND)
  if(RUNNING_CGAL_AUTO_TEST)
    # Just to avoid a warning from CMake if that variable is set on the command line...
  endif()

  message(STATUS "NOTICE: This demo needs Qt5, and will not be compiled.")
endif(
  CGAL_Qt5_FOUND
  AND CGAL_ImageIO_FOUND)
