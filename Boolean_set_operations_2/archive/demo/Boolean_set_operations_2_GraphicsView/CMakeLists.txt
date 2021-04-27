# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.1...3.20)
project(Boolean_set_operations_2_GraphicsView_Demo)

if(NOT POLICY CMP0070 AND POLICY CMP0053)
  # Only set CMP0053 to OLD with CMake<3.10, otherwise there is a warning.
  cmake_policy(SET CMP0053 OLD)
endif()

find_package(CGAL COMPONENTS Qt5 Core)

set(QT_USE_QTMAIN TRUE)

find_package(Qt5 QUIET COMPONENTS Widgets Script Svg)

include_directories(BEFORE ${CMAKE_CURRENT_SOURCE_DIR}/include)

if(CGAL_Qt5_FOUND
   AND Qt5_FOUND
   AND CGAL_Core_FOUND)
  set(CMAKE_AUTOMOC ON)
  # UI files (Qt Designer files)
  qt5_wrap_ui(CDT_UI_FILES boolean_operations_2.ui)

  # qrc files (resources files, that contain icons, at least)
  qt5_add_resources(CGAL_Qt5_RESOURCE_FILES
                    ${CMAKE_CURRENT_SOURCE_DIR}/boolean_operations_2.qrc)

  # use the Qt MOC preprocessor on classes that derives from QObject

  # The executable itself.
  add_executable(
    boolean_operations_2
    ${CMAKE_CURRENT_SOURCE_DIR}/boolean_operations_2.cpp ${CGAL_Qt5_MOC_FILES}
    ${CDT_UI_FILES} ${CGAL_Qt5_RESOURCE_FILES})
  # Link with Qt libraries
  target_link_libraries(boolean_operations_2 Qt5::Widgets Qt5::Script Qt5::Svg)

  # Link with CGAL
  target_link_libraries(boolean_operations_2 ${CGAL_LIBRARIES}
                        ${CGAL_3RD_PARTY_LIBRARIES})

else()

  message(
    STATUS "NOTICE: This demo requires CGAL and Qt5, and will not be compiled.")

endif()
