message("Configuring libCGAL_Qt5")

if(POLICY CMP0043)
  cmake_policy(SET CMP0043 OLD)
endif()

find_package(Qt5 QUIET COMPONENTS OpenGL Svg QUIET)
find_package(OpenGL QUIET)

set(CGAL_Qt5_MISSING_DEPS "")
if(NOT Qt5OpenGL_FOUND)
  set(CGAL_Qt5_MISSING_DEPS "Qt5OpenGL")
endif()
if(NOT Qt5Svg_FOUND)
  set(CGAL_Qt5_MISSING_DEPS "${CGAL_Qt5_MISSING_DEPS} Qt5Svg")
endif()
if(NOT Qt5_FOUND)
  set(CGAL_Qt5_MISSING_DEPS "${CGAL_Qt5_MISSING_DEPS} Qt5")
endif()
if(NOT OPENGL_FOUND)
  set(CGAL_Qt5_MISSING_DEPS "${CGAL_Qt5_MISSING_DEPS} OpenGL")
endif()

if(CGAL_Qt5_MISSING_DEPS)
  message(STATUS "libCGAL_Qt5 is missing the dependencies: ${CGAL_Qt5_MISSING_DEPS} cannot be configured.")
  return()
endif()

include_directories (BEFORE ../../include)

message( STATUS "USING Qt5_VERSION = '${Qt5Core_VERSION_STRING}'" )
if(COMMAND add_config_flag)
  set( CGAL_HAS_QT5 TRUE )
  add_config_flag( CGAL_HAS_QT5 )
endif()

get_property(QT_UIC_EXECUTABLE TARGET Qt5::uic PROPERTY LOCATION)
message( STATUS "OpenGL include:      ${OPENGL_INCLUDE_DIR}" )
message( STATUS "OpenGL libraries:    ${OPENGL_LIBRARIES}" )
message( STATUS "OpenGL definitions:  ${OPENGL_DEFINITIONS}" )
message( STATUS "Qt5Core include:     ${Qt5Core_INCLUDE_DIRS}" )
message( STATUS "Qt5 libraries:       ${Qt5Core_LIBRARIES} ${Qt5Gui_LIBRARIES} ${Qt5Svg_LIBRARIES} ${Qt5OpenGL_LIBRARIES}" )
message( STATUS "Qt5Core definitions: ${Qt5Core_DEFINITIONS}" )
message( STATUS "moc executable:      ${QT_MOC_EXECUTABLE}" )
message( STATUS "uic executable:      ${QT_UIC_EXECUTABLE}" )

cache_set(CGAL_Qt5_3RD_PARTY_INCLUDE_DIRS ${QT_INCLUDE_DIR} ${OPENGL_INCLUDE_DIR} )
cache_set(CGAL_Qt5_3RD_PARTY_LIBRARIES    ${QT_LIBRARIES}   ${OPENGL_LIBRARIES}   )
cache_set(CGAL_Qt5_3RD_PARTY_DEFINITIONS  ${QT_DEFINITIONS} ${OPENGL_DEFINITIONS} )

set(mocfiles "")
set(RESOURCE_FILES "")

foreach (package ${CGAL_CONFIGURED_PACKAGES} )
  file(GLOB PACKAGE_QTMOC_FILES "${package}/src/CGAL_Qt5/*.qtmoc.cmake")
  foreach(package_qtmoc_file ${PACKAGE_QTMOC_FILES})
    # includes 'moccing' for sources/headers in package + collects lists of moc-files for dependency (to properly build the lib)
    include(${package_qtmoc_file}) 
    #        message(STATUS QTMOC------------FILE: ${package_qtmoc_file})      
  endforeach()
endforeach()

foreach(mocfile ${mocfiles})
  list(APPEND additional_files ${mocfile})
endforeach()
foreach(resfile ${RESOURCE_FILES})
  list(APPEND additional_files ${resfile})
endforeach()

#    message(STATUS "Additional input files: ${additional_files}")

use_essential_libs()

include_directories( SYSTEM ${CGAL_3RD_PARTY_INCLUDE_DIRS} ${CGAL_Qt5_3RD_PARTY_INCLUDE_DIRS} )

link_directories    ( ${CGAL_LIBRARIES_DIR} ${CGAL_3RD_PARTY_LIBRARIES_DIRS} )

collect_cgal_library( CGAL_Qt5 "${additional_files}")

qt5_use_modules(CGAL_Qt5 OpenGL Svg)

add_dependencies( CGAL_Qt5 CGAL )

# CGAL_Qt5 only depends on CGAL in DEBUG, but we still link it in
# both build types.
target_link_libraries( CGAL_Qt5 CGAL ${CGAL_3RD_PARTY_LIBRARIES} ${CGAL_Qt5_3RD_PARTY_LIBRARIES} )

add_definitions ( ${CGAL_3RD_PARTY_DEFINITIONS} ${CGAL_Qt5_3RD_PARTY_DEFINITIONS} )

if($ENV{CGAL_FAKE_PUBLIC_RELEASE})
  add_definitions( -DCGAL_FAKE_PUBLIC_RELEASE )
endif()

message("libCGAL_Qt5 is configured")
