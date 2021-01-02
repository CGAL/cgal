# This is the CMake script for compiling the CGALRemeshing.

project( CGALRemeshing )

# Find includes in corresponding build directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)
# Instruct CMake to run moc automatically when needed.
set(CMAKE_AUTOMOC ON)
cmake_minimum_required(VERSION 3.1)
if(NOT POLICY CMP0070 AND POLICY CMP0053)
  # Only set CMP0053 to OLD with CMake<3.10, otherwise there is a warning.
  cmake_policy(SET CMP0053 OLD)
endif()

if(POLICY CMP0071)
  cmake_policy(SET CMP0071 NEW)
endif()

# Include this package's headers first
include_directories( BEFORE ./ ./include ../../include )

# Find CGAL and CGAL Qt5
find_package(CGAL COMPONENTS Qt5)

# Find Qt5 itself
find_package(Qt5 QUIET COMPONENTS Xml Script OpenGL Gui Svg)

include( ${CGAL_USE_FILE} )


if(CGAL_FOUND AND CGAL_Qt5_FOUND AND Qt5_FOUND)

  qt5_wrap_ui( UI_FILES Main_window.ui Isotropic_parameters.ui MinAngle_parameters.ui)

  include(AddFileDependencies)

  qt5_generate_moc( "Main_window.h" "${CMAKE_CURRENT_BINARY_DIR}/Main_window_moc.cpp" )
  add_file_dependencies( Main_window_moc.cpp "${CMAKE_CURRENT_SOURCE_DIR}/Main_window.h" )

  qt5_generate_moc( "Viewer.h" "${CMAKE_CURRENT_BINARY_DIR}/Viewer_moc.cpp" )
  add_file_dependencies( Viewer_moc.cpp "${CMAKE_CURRENT_SOURCE_DIR}/Viewer.h" )

  qt5_generate_moc( "Scene.h" "${CMAKE_CURRENT_BINARY_DIR}/Scene_moc.cpp" )
  add_file_dependencies( Scene_moc.cpp "${CMAKE_CURRENT_SOURCE_DIR}/Scene.h" )

  qt5_generate_moc( "Isotropic_parameters.h" "${CMAKE_CURRENT_BINARY_DIR}/Isotropic_parameters_moc.cpp" )
  add_file_dependencies( Isotropic_parameters_moc.cpp "${CMAKE_CURRENT_SOURCE_DIR}/Isotropic_parameters.h" )

  qt5_generate_moc( "MinAngle_parameters.h" "${CMAKE_CURRENT_BINARY_DIR}/MinAngle_parameters_moc.cpp" )
  add_file_dependencies( MinAngle_parameters_moc.cpp "${CMAKE_CURRENT_SOURCE_DIR}/MinAngle_parameters.h" )

  qt5_add_resources ( CGAL_Qt5_RESOURCE_FILES CGAL_remeshing.qrc )

  add_file_dependencies( CGAL_remeshing.cpp "${CMAKE_CURRENT_BINARY_DIR}/Main_window_moc.cpp"
                                       "${CMAKE_CURRENT_BINARY_DIR}/Viewer_moc.cpp"
                                       "${CMAKE_CURRENT_BINARY_DIR}/Scene_moc.cpp" 
	                       "${CMAKE_CURRENT_BINARY_DIR}/Isotropic_parameters_moc.cpp"
		       "${CMAKE_CURRENT_BINARY_DIR}/MinAngle_parameters_moc.cpp")

  add_executable  ( CGALRemeshing CGAL_remeshing.cpp ${UI_FILES} ${CGAL_Qt5_RESOURCE_FILES} 
    #${CGAL_Qt5_MOC_FILES}
    )
  # Link with Qt libraries
  target_link_libraries( CGALRemeshing PRIVATE
    Qt5::OpenGL Qt5::Gui Qt5::Xml
    CGAL::CGAL
    CGAL::CGAL_Qt5
    )

  add_to_cached_list( CGAL_EXECUTABLE_TARGETS CGALRemeshing )

  include(${CGAL_MODULES_DIR}/CGAL_add_test.cmake)
  cgal_add_compilation_test(CGALRemeshing)

else (CGAL_FOUND AND CGAL_Qt5_FOUND AND Qt5_FOUND)

  set(RR_MISSING_DEPS "")

  if(NOT CGAL_FOUND)
    set(RR_MISSING_DEPS "the CGAL library, ${RR_MISSING_DEPS}")
  endif()

  if(NOT CGAL_Qt5_FOUND)
    set(RR_MISSING_DEPS "CGAL_Qt5, ${RR_MISSING_DEPS}")
  endif()

  if(NOT Qt5_FOUND)
    set(RR_MISSING_DEPS "Qt5, ${RR_MISSING_DEPS}")
  endif()

  message(STATUS "NOTICE: This demo requires ${RR_MISSING_DEPS}and will not be compiled.")

endif (CGAL_FOUND AND CGAL_Qt5_FOUND AND Qt5_FOUND )
