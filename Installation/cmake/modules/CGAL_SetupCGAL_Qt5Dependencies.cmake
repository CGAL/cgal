#.rst:
# CGAL_SetupCGAL_Qt5Dependencies
# ------------------------------
#
# The module searchs for the dependencies of the `CGAL_Qt5` library:
#   - the `Qt5` libraries
#
# by calling
#
# .. code-block:: cmake
#
#    find_package(Qt5 QUIET COMPONENTS OpenGL Svg)
#
# and defines the variable :variable:`CGAL_Qt5_FOUND` and the function
# :command:`CGAL_setup_CGAL_Qt5_dependencies`.
#

if(CGAL_SetupCGAL_Qt5Dependencies_included)
  return()
endif()
set(CGAL_SetupCGAL_Qt5Dependencies_included TRUE)

#.rst:
# Used Modules
# ^^^^^^^^^^^^
#   - :module:`Qt5Config`
find_package(Qt5 QUIET COMPONENTS OpenGL Svg Xml)

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
if(NOT EXISTS ${CGAL_GRAPHICSVIEW_PACKAGE_DIR}/include/CGAL/Qt/GraphicsItem.h)
  set(CGAL_Qt5_MISSING_DEPS "${CGAL_Qt5_MISSING_DEPS} <CGAL/Qt/*.h> headers")
endif()

#.rst:
# Result Variables
# ^^^^^^^^^^^^^^^^
#
# .. variable:: CGAL_Qt5_FOUND
#
#    Set to `TRUE` if the dependencies of `CGAL_Qt5` were found.
#
if(NOT CGAL_Qt5_MISSING_DEPS)
  set(CGAL_Qt5_FOUND TRUE)
  set_property(GLOBAL PROPERTY CGAL_Qt5_FOUND TRUE)

  include(${CMAKE_CURRENT_LIST_DIR}/CGAL_Qt5_moc_and_resource_files.cmake)

  if(CGAL_HEADER_ONLY AND (WITH_demos OR WITH_examples OR WITH_tests OR NOT CGAL_BUILDING_LIBS) AND NOT TARGET CGAL_Qt5_moc_and_resources)
    add_library(CGAL_Qt5_moc_and_resources STATIC 
      ${_CGAL_Qt5_MOC_FILES_private} 
        ${CGAL_GRAPHICSVIEW_PACKAGE_DIR}/include/CGAL/Qt/GraphicsViewNavigation.h
        ${CGAL_GRAPHICSVIEW_PACKAGE_DIR}/include/CGAL/Qt/DemosMainWindow.h
        ${CGAL_GRAPHICSVIEW_PACKAGE_DIR}/include/CGAL/Qt/GraphicsItem.h
        ${CGAL_GRAPHICSVIEW_PACKAGE_DIR}/include/CGAL/Qt/GraphicsViewInput.h
        ${CGAL_GRAPHICSVIEW_PACKAGE_DIR}/include/CGAL/Qt/camera.h
        ${CGAL_GRAPHICSVIEW_PACKAGE_DIR}/include/CGAL/Qt/frame.h
        ${CGAL_GRAPHICSVIEW_PACKAGE_DIR}/include/CGAL/Qt/keyFrameInterpolator.h
        ${CGAL_GRAPHICSVIEW_PACKAGE_DIR}/include/CGAL/Qt/manipulatedCameraFrame.h
        ${CGAL_GRAPHICSVIEW_PACKAGE_DIR}/include/CGAL/Qt/manipulatedFrame.h
        ${CGAL_GRAPHICSVIEW_PACKAGE_DIR}/include/CGAL/Qt/qglviewer.h
        ${CGAL_GRAPHICSVIEW_PACKAGE_DIR}/include/CGAL/Qt/image_interface.h
        ${_CGAL_Qt5_UI_FILES}
      ${_CGAL_Qt5_RESOURCE_FILES_private})
    target_include_directories( CGAL_Qt5_moc_and_resources PUBLIC ${CMAKE_CURRENT_BINARY_DIR})
    set_target_properties(CGAL_Qt5_moc_and_resources PROPERTIES
      POSITION_INDEPENDENT_CODE TRUE
      EXCLUDE_FROM_ALL TRUE
      AUTOMOC TRUE)
    target_link_libraries(CGAL_Qt5_moc_and_resources CGAL::CGAL Qt5::Widgets Qt5::OpenGL Qt5::Svg Qt5::Xml)

    add_library(CGAL::CGAL_Qt5_moc_and_resources ALIAS CGAL_Qt5_moc_and_resources)
    add_library(CGAL::Qt5_moc_and_resources ALIAS CGAL_Qt5_moc_and_resources)
  endif()

endif()

#get_property(QT_UIC_EXECUTABLE TARGET Qt5::uic PROPERTY LOCATION)
#message( STATUS "Qt5Core include:     ${Qt5Core_INCLUDE_DIRS}" )
#message( STATUS "Qt5 libraries:       ${Qt5Core_LIBRARIES} ${Qt5Gui_LIBRARIES} ${Qt5Svg_LIBRARIES} ${Qt5OpenGL_LIBRARIES}" )
#message( STATUS "Qt5Core definitions: ${Qt5Core_DEFINITIONS}" )
#message( STATUS "moc executable:      ${QT_MOC_EXECUTABLE}" )
#message( STATUS "uic executable:      ${QT_UIC_EXECUTABLE}" )

#.rst:
#
# Provided Functions
# ^^^^^^^^^^^^^^^^^^
#
# .. command:: CGAL_setup_CGAL_Qt5_dependencies
#
#   Link the target with the dependencies of `CGAL_Qt5`::
#
#     CGAL_setup_CGAL_Qt5_dependencies( target [INTERFACE] )
#
#   If the option ``INTERFACE`` is passed, the dependencies are
#   added using :command:`target_link_libraries` with the ``INTERFACE``
#   keyword, or ``PUBLIC`` otherwise.
#
function(CGAL_setup_CGAL_Qt5_dependencies target)
  if(ARGV1 STREQUAL INTERFACE)
    set(keyword INTERFACE)
  else()
    set(keyword PUBLIC)
  endif()

  if($ENV{CGAL_FAKE_PUBLIC_RELEASE})
    target_compile_definitions( ${target} ${keyword} CGAL_FAKE_PUBLIC_RELEASE=1 )
  endif()
  target_link_libraries( ${target} ${keyword} CGAL::CGAL)
  if(CGAL_HEADER_ONLY)
    target_link_libraries( ${target} ${keyword} CGAL::Qt5_moc_and_resources)
  endif()
  target_link_libraries( ${target} ${keyword} Qt5::OpenGL Qt5::Svg Qt5::Xml)

  # Remove -Wdeprecated-copy, for g++ >= 9.0, because Qt5, as of
  # version 5.12, has a lot of [-Wdeprecated-copy] warnings.
  if( CMAKE_CXX_COMPILER_ID STREQUAL "GNU"
      AND NOT CMAKE_CXX_COMPILER_VERSION VERSION_LESS "9" )
    target_compile_options( ${target} ${keyword} "-Wno-deprecated-copy" "-Wno-cast-function-type" )
  endif()

endfunction()

