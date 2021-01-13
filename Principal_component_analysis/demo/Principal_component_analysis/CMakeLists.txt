# This is the CMake script for compiling the PCA demo.

cmake_minimum_required(VERSION 3.1...3.15)
project(Principal_component_analysis_Demo)

if(NOT POLICY CMP0070 AND POLICY CMP0053)
  # Only set CMP0053 to OLD with CMake<3.10, otherwise there is a warning.
  cmake_policy(SET CMP0053 OLD)
endif()

if(POLICY CMP0071)
  cmake_policy(SET CMP0071 NEW)
endif()

include_directories(./)

# Find CGAL and CGAL Qt5
find_package(CGAL REQUIRED OPTIONAL_COMPONENTS Qt5)

find_package(Eigen3 3.1.0) #(requires 3.1.0 or greater)
include(CGAL_Eigen_support)
if(NOT TARGET CGAL::Eigen_support)
  message(
    STATUS
      "NOTICE: This project requires the Eigen library, and will not be compiled."
  )
  return()
endif()

# Find Qt5 itself
find_package(Qt5 QUIET COMPONENTS Xml Script OpenGL)

if(CGAL_Qt5_FOUND AND Qt5_FOUND)

  set(CMAKE_INCLUDE_CURRENT_DIR ON)
  qt5_wrap_ui(UI_FILES MainWindow.ui)

  include(AddFileDependencies)

  qt5_generate_moc("MainWindow.h"
                   "${CMAKE_CURRENT_BINARY_DIR}/MainWindow_moc.cpp")
  add_file_dependencies(MainWindow_moc.cpp
                        "${CMAKE_CURRENT_SOURCE_DIR}/MainWindow.h")

  qt5_generate_moc("Viewer.h" "${CMAKE_CURRENT_BINARY_DIR}/Viewer_moc.cpp")
  add_file_dependencies(Viewer_moc.cpp "${CMAKE_CURRENT_SOURCE_DIR}/Viewer.h")

  qt5_add_resources(CGAL_Qt5_RESOURCE_FILES PCA_demo.qrc)

  add_file_dependencies(
    PCA_demo.cpp "${CMAKE_CURRENT_BINARY_DIR}/MainWindow_moc.cpp"
    "${CMAKE_CURRENT_BINARY_DIR}/Viewer_moc.cpp")
  add_executable(PCA_demo PCA_demo.cpp ${UI_FILES} ${CGAL_Qt5_RESOURCE_FILES}
                          ${CGAL_Qt5_MOC_FILES})

  target_link_libraries(PCA_demo PRIVATE CGAL::CGAL CGAL::CGAL_Qt5
                                         CGAL::Eigen_support Qt5::Gui)

  add_to_cached_list(CGAL_EXECUTABLE_TARGETS PCA_demo)

  include(${CGAL_MODULES_DIR}/CGAL_add_test.cmake)
  cgal_add_compilation_test(PCA_demo)

else(CGAL_Qt5_FOUND AND Qt5_FOUND)

  set(PCA_MISSING_DEPS "")

  if(NOT CGAL_Qt5_FOUND)
    set(PCA_MISSING_DEPS "the CGAL Qt5 library, ${PCA_MISSING_DEPS}")
  endif()

  if(NOT Qt5_FOUND)
    set(PCA_MISSING_DEPS "Qt5, ${PCA_MISSING_DEPS}")
  endif()

  message(
    STATUS
      "NOTICE: This demo requires ${PCA_MISSING_DEPS} and will not be compiled."
  )

endif(CGAL_Qt5_FOUND AND Qt5_FOUND)
