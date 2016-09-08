message("Configuring libCGAL_Qt5")

include(CGAL_SetupCGAL_Qt5Dependencies)

if(CGAL_Qt5_MISSING_DEPS)
  message(STATUS "libCGAL_Qt5 is missing the dependencies: ${CGAL_Qt5_MISSING_DEPS} cannot be configured.")
  return()
endif()

message( STATUS "USING Qt5_VERSION = '${Qt5Core_VERSION_STRING}'" )

include(Use_CGAL_Qt5_headers)

list(APPEND additional_files ${CGAL_Qt5_MOC_FILES} ${CGAL_Qt5_RESOURCE_FILES})

collect_cgal_library( CGAL_Qt5 "${additional_files}")

if(CGAL_HEADER_ONLY)
  set(keyword "INTERFACE")
endif()

CGAL_setup_CGAL_Qt5_dependencies( CGAL_Qt5 ${keyword})

if (CGAL_HEADER_ONLY)
  install(DIRECTORY "../../include/CGAL/Qt/" DESTINATION "${CGAL_INSTALL_INC_DIR}/CGAL/include/CGAL/Qt")
  install(DIRECTORY "../../demo/resources/" DESTINATION "${CGAL_INSTALL_INC_DIR}/CGAL/demo/resources")
  install(DIRECTORY "../../demo/icons/" DESTINATION "${CGAL_INSTALL_INC_DIR}/CGAL/demo/icons")

endif()

message("libCGAL_Qt5 is configured")
