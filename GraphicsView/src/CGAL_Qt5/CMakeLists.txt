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

collect_cgal_library( CGAL_Qt5 "")

set(keyword "INTERFACE")

CGAL_setup_CGAL_Qt5_dependencies( CGAL_Qt5 ${keyword})

if(COMMAND add_config_flag)
  set( CGAL_HAS_QT5 TRUE )
  add_config_flag( CGAL_HAS_QT5 )
endif()

install(DIRECTORY "../../include/CGAL/Qt/" DESTINATION "${CGAL_INSTALL_INC_DIR}/CGAL/Qt" COMPONENT CGAL_Qt5)
install(DIRECTORY "../../demo/resources/" DESTINATION "${CGAL_INSTALL_CMAKE_DIR}/demo/resources" COMPONENT CGAL_Qt5)
install(DIRECTORY "../../demo/icons/" DESTINATION "${CGAL_INSTALL_CMAKE_DIR}/demo/icons" COMPONENT CGAL_Qt5)

message("libCGAL_Qt5 is configured")
