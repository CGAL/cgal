# This is the CMake script for compiling the PCA demo.

project( Principal_component_analysis_Demo )

cmake_minimum_required(VERSION 3.1)
if(NOT POLICY CMP0070 AND POLICY CMP0053)
  # Only set CMP0053 to OLD with CMake<3.10, otherwise there is a warning.
  cmake_policy(SET CMP0053 OLD)
endif()

foreach(INCDIR ../../include ../../../STL_Extension/include ../../../GraphicsView/include ../../../filtered_kernel/include )
  if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${INCDIR}")
    include_directories (BEFORE "${CMAKE_CURRENT_SOURCE_DIR}/${INCDIR}")
  endif()
endforeach()

include_directories( ./ )

# Find CGAL and CGAL Qt5
find_package(CGAL COMPONENTS Qt5)
include( ${CGAL_USE_FILE} )

find_package(Eigen3 3.1.0) #(requires 3.1.0 or greater)
if (EIGEN3_FOUND)
  include( ${EIGEN3_USE_FILE} )
else()
  message(STATUS "NOTICE: This project requires the Eigen library, and will not be compiled.")
  return()
endif()

# Find Qt5 itself
find_package(Qt5 QUIET COMPONENTS Xml Script OpenGL)

if(CGAL_Qt5_FOUND AND Qt5_FOUND )

  qt5_wrap_ui( UI_FILES MainWindow.ui )

  include(AddFileDependencies)

  qt5_generate_moc( "MainWindow.h" "${CMAKE_CURRENT_BINARY_DIR}/MainWindow_moc.cpp" )
  add_file_dependencies( MainWindow_moc.cpp "${CMAKE_CURRENT_SOURCE_DIR}/MainWindow.h" )

  qt5_generate_moc( "Viewer.h" "${CMAKE_CURRENT_BINARY_DIR}/Viewer_moc.cpp" )
  add_file_dependencies( Viewer_moc.cpp "${CMAKE_CURRENT_SOURCE_DIR}/Viewer.h" )

  qt5_add_resources ( CGAL_Qt5_RESOURCE_FILES PCA_demo.qrc )

  add_file_dependencies( PCA_demo.cpp "${CMAKE_CURRENT_BINARY_DIR}/MainWindow_moc.cpp"
    "${CMAKE_CURRENT_BINARY_DIR}/Viewer_moc.cpp" )
  add_executable  ( PCA_demo PCA_demo.cpp ${UI_FILES} ${CGAL_Qt5_RESOURCE_FILES} ${CGAL_Qt5_MOC_FILES})

  target_link_libraries( PCA_demo PRIVATE
    CGAL::CGAL CGAL::CGAL_Qt5 Qt5::Gui)

  add_to_cached_list( CGAL_EXECUTABLE_TARGETS PCA_demo )

  include(${CGAL_MODULES_DIR}/CGAL_add_test.cmake)
  cgal_add_compilation_test(PCA_demo)

else (CGAL_Qt5_FOUND AND Qt5_FOUND )

  set(PCA_MISSING_DEPS "")

  if(NOT CGAL_Qt5_FOUND)
    set(PCA_MISSING_DEPS "the CGAL Qt5 library, ${PCA_MISSING_DEPS}")
  endif()

  if(NOT Qt5_FOUND)
    set(PCA_MISSING_DEPS "Qt5, ${PCA_MISSING_DEPS}")
  endif()

  
  message(STATUS "NOTICE: This demo requires ${PCA_MISSING_DEPS} and will not be compiled.")

endif (CGAL_Qt5_FOUND AND Qt5_FOUND )
