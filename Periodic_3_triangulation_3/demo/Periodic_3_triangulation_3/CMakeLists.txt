# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.1...3.15)
project(Periodic_3_triangulation_3_Demo)

if(NOT POLICY CMP0070 AND POLICY CMP0053)
  # Only set CMP0053 to OLD with CMake<3.10, otherwise there is a warning.
  cmake_policy(SET CMP0053 OLD)
endif()

if(POLICY CMP0071)
  cmake_policy(SET CMP0071 NEW)
endif()

# Find CGAL
find_package(CGAL REQUIRED OPTIONAL_COMPONENTS Qt5)

# Find Qt5 itself
find_package(Qt5 QUIET COMPONENTS Xml OpenGL Help Core)

if(Qt5_FOUND)
  add_definitions(-DQT_NO_KEYWORDS)
  set(CMAKE_INCLUDE_CURRENT_DIR ON)
endif(Qt5_FOUND)

if(Qt5Help_VERSION VERSION_LESS 5.12)
  set(CGAL_QCOLLECTIONGENERATOR_TARGET Qt5::qcollectiongenerator)
else()
  set(CGAL_QCOLLECTIONGENERATOR_TARGET Qt5::qhelpgenerator)
endif()

if(CGAL_Qt5_FOUND
   AND Qt5_FOUND
   AND TARGET ${CGAL_QCOLLECTIONGENERATOR_TARGET})

  # UI files (Qt Designer files)
  qt5_wrap_ui(UI_FILES MainWindow.ui)

  # qrc files (resource files)
  qt5_add_resources(CGAL_Qt5_RESOURCE_FILES ./Periodic_3_triangulation_3.qrc)

  # use the Qt MOC preprocessor on classes that derive from QObject
  qt5_generate_moc("Scene.h" "${CMAKE_CURRENT_BINARY_DIR}/moc_Scene.cpp")
  qt5_generate_moc("MainWindow.h"
                   "${CMAKE_CURRENT_BINARY_DIR}/moc_MainWindow.cpp")
  qt5_generate_moc("Viewer.h" "${CMAKE_CURRENT_BINARY_DIR}/moc_Viewer.cpp")

  # generate QtAssistant collection file
  add_custom_command(
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/Periodic_3_triangulation_3.qhc
    DEPENDS Periodic_3_triangulation_3.qhp Periodic_3_triangulation_3.qhcp
    COMMAND
      ${CGAL_QCOLLECTIONGENERATOR_TARGET}
      ${CMAKE_CURRENT_SOURCE_DIR}/Periodic_3_triangulation_3.qhcp -o
      ${CMAKE_CURRENT_BINARY_DIR}/Periodic_3_triangulation_3.qhc)

  # The executable itself
  add_executable(
    periodic_3_triangulation_3_demo
    Scene.cpp
    moc_Scene.cpp
    Viewer.cpp
    moc_Viewer.cpp
    periodic_3_triangulation_3_demo.cpp
    MainWindow.ui
    moc_MainWindow.cpp
    ${UI_FILES}
    ${CGAL_Qt5_RESOURCE_FILES}
    ${CGAL_Qt5_MOC_FILES}
    Periodic_3_triangulation_3.qhc)

  add_to_cached_list(CGAL_EXECUTABLE_TARGETS periodic_3_triangulation_3_demo)

  target_link_libraries(periodic_3_triangulation_3_demo
                        PRIVATE CGAL::CGAL CGAL::CGAL_Qt5 Qt5::OpenGL)

  include(${CGAL_MODULES_DIR}/CGAL_add_test.cmake)
  cgal_add_compilation_test(periodic_3_triangulation_3_demo)
else()

  set(PERIODIC_TRIANGULATION_MISSING_DEPS "")

  if(NOT CGAL_Qt5_FOUND)
    set(PERIODIC_TRIANGULATION_MISSING_DEPS
        "the CGAL Qt5 library, ${PERIODIC_TRIANGULATION_MISSING_DEPS}")
  endif()

  if(NOT Qt5_FOUND)
    set(PERIODIC_TRIANGULATION_MISSING_DEPS
        "Qt5, ${PERIODIC_TRIANGULATION_MISSING_DEPS}")
  endif()

  if(NOT TARGET ${CGAL_QCOLLECTIONGENERATOR_TARGET})
    set(PERIODIC_TRIANGULATION_MISSING_DEPS
        "${CGAL_QCOLLECTIONGENERATOR_TARGET}, ${PERIODIC_TRIANGULATION_MISSING_DEPS}"
    )
  endif()

  message(
    STATUS
      "NOTICE: This demo requires ${PERIODIC_TRIANGULATION_MISSING_DEPS}and will not be compiled."
  )

endif()
