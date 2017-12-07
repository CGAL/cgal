message("Configuring libCGAL_Qt5")

if(NOT POLICY CMP0070 AND POLICY CMP0053)
  # Only set CMP0053 to OLD with CMake<3.10, otherwise there is a warning.
  cmake_policy(SET CMP0053 OLD)
endif()

include(CGAL_SetupCGAL_Qt5Dependencies)

if(CGAL_Qt5_MISSING_DEPS)
  message(STATUS "libCGAL_Qt5 is missing the dependencies: ${CGAL_Qt5_MISSING_DEPS} cannot be configured.")
  return()
endif()

message( STATUS "USING Qt5_VERSION = '${Qt5Core_VERSION_STRING}'" )

if(NOT CGAL_HEADER_ONLY)
  collect_cgal_library( CGAL_Qt5 "${_CGAL_Qt5_MOC_FILES_private};${_CGAL_Qt5_RESOURCE_FILES_private}")
else()
  collect_cgal_library( CGAL_Qt5 "")
endif()

if(CGAL_HEADER_ONLY)
  set(keyword "INTERFACE")
endif()

CGAL_setup_CGAL_Qt5_dependencies( CGAL_Qt5 ${keyword})

if(COMMAND add_config_flag)
  set( CGAL_HAS_QT5 TRUE )
  add_config_flag( CGAL_HAS_QT5 )
endif()

if (CGAL_HEADER_ONLY)
  install(DIRECTORY "../../include/CGAL/Qt/" DESTINATION "${CGAL_INSTALL_INC_DIR}/CGAL/include/CGAL/Qt" COMPONENT CGAL_Qt5)
  install(DIRECTORY "../../demo/resources/" DESTINATION "${CGAL_INSTALL_CMAKE_DIR}/demo/resources" COMPONENT CGAL_Qt5)
  install(DIRECTORY "../../demo/icons/" DESTINATION "${CGAL_INSTALL_CMAKE_DIR}/demo/icons" COMPONENT CGAL_Qt5)
endif()

message("libCGAL_Qt5 is configured")
