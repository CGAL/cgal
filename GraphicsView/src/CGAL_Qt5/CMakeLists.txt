message("Configuring libCGAL_Qt5")

if(NOT POLICY CMP0070 AND POLICY CMP0053)
  # Only set CMP0053 to OLD with CMake<3.10, otherwise there is a warning.
  cmake_policy(SET CMP0053 OLD)
endif()

include(CGAL_SetupCGAL_Qt5Dependencies)

if(CGAL_Qt5_MISSING_DEPS)
  message(
    STATUS
      "libCGAL_Qt5 is missing the dependencies: ${CGAL_Qt5_MISSING_DEPS} cannot be configured."
  )
  return()
endif()

message(STATUS "USING Qt5_VERSION = '${Qt5Core_VERSION_STRING}'")
if(NOT CGAL_HEADER_ONLY)
  collect_cgal_library(
    CGAL_Qt5
    "${_CGAL_Qt5_MOC_FILES_private};${_CGAL_Qt5_RESOURCE_FILES_private};${_CGAL_Qt5_UI_FILES}"
  )
  target_include_directories(CGAL_Qt5 PRIVATE ${CMAKE_CURRENT_BINARY_DIR})
else()
  collect_cgal_library(CGAL_Qt5 "")
endif()

if(CGAL_HEADER_ONLY)
  set(keyword "INTERFACE")
endif()

cgal_setup_cgal_qt5_dependencies(CGAL_Qt5 ${keyword})

if(COMMAND add_config_flag)
  set(CGAL_HAS_QT5 TRUE)
  add_config_flag(CGAL_HAS_QT5)
endif()

install(
  DIRECTORY "../../include/CGAL/Qt/"
  DESTINATION "${CGAL_INSTALL_INC_DIR}/CGAL/Qt"
  COMPONENT CGAL_Qt5)
if(CGAL_HEADER_ONLY)
  if(EXISTS ${CGAL_MODULES_DIR}/demo/resources)
    install(
      DIRECTORY "${CGAL_MODULES_DIR}/demo/resources/"
      DESTINATION "${CGAL_INSTALL_CMAKE_DIR}/demo/resources"
      COMPONENT CGAL_Qt5)
    install(
      DIRECTORY "${CGAL_MODULES_DIR}/demo/icons/"
      DESTINATION "${CGAL_INSTALL_CMAKE_DIR}/demo/icons"
      COMPONENT CGAL_Qt5)
  else()
    install(
      DIRECTORY "../../demo/resources/"
      DESTINATION "${CGAL_INSTALL_CMAKE_DIR}/demo/resources"
      COMPONENT CGAL_Qt5)
    install(
      DIRECTORY "../../demo/icons/"
      DESTINATION "${CGAL_INSTALL_CMAKE_DIR}/demo/icons"
      COMPONENT CGAL_Qt5)
  endif()
endif()

message("libCGAL_Qt5 is configured")
