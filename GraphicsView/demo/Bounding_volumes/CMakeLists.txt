# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.1...3.15)
project(Bounding_volumes_Demo)

if(NOT POLICY CMP0070 AND POLICY CMP0053)
  # Only set CMP0053 to OLD with CMake<3.10, otherwise there is a warning.
  cmake_policy(SET CMP0053 OLD)
endif()

if(POLICY CMP0071)
  cmake_policy(SET CMP0071 NEW)
endif()

find_package(CGAL REQUIRED OPTIONAL_COMPONENTS Qt5)

find_package(Qt5 QUIET COMPONENTS Xml Script OpenGL Svg)
include_directories(BEFORE ./include)

if(CGAL_Qt5_FOUND AND Qt5_FOUND)

  set(CMAKE_INCLUDE_CURRENT_DIR ON)
  add_definitions(-DQT_NO_KEYWORDS)
  set(CMAKE_AUTOMOC ON)

  #----------------------------------------------
  # The "Bounding volumes" demo: Bounding_volumes
  #----------------------------------------------
  # UI files (Qt Designer files)
  qt5_wrap_ui(DT_UI_FILES Bounding_volumes.ui)

  # qrc files (resources files, that contain icons, at least)
  qt5_add_resources(CGAL_Qt5_RESOURCE_FILES ./Bounding_volumes.qrc)

  # use the Qt MOC preprocessor on classes that derives from QObject

  # The executable itself.
  add_executable(
    Bounding_volumes Bounding_volumes.cpp ${DT_UI_FILES}
                     ${CGAL_Qt5_RESOURCE_FILES} ${CGAL_Qt5_MOC_FILES})

  add_to_cached_list(CGAL_EXECUTABLE_TARGETS Bounding_volumes)

  target_link_libraries(Bounding_volumes PRIVATE CGAL::CGAL CGAL::CGAL_Qt5
                                                 Qt5::Gui)

  include(${CGAL_MODULES_DIR}/CGAL_add_test.cmake)
  cgal_add_compilation_test(Bounding_volumes)

  include(${CGAL_MODULES_DIR}/CGAL_add_test.cmake)
  cgal_add_compilation_test(Bounding_volumes)

else()

  message(
    STATUS "NOTICE: This demo requires CGAL and Qt5, and will not be compiled.")

endif()
