#.rst:
# CGAL_SetupCGAL_Qt6Dependencies
# ------------------------------
#
# The module searches for the dependencies of the `CGAL_Qt6` library:
#   - the `Qt6` libraries
#
# by calling
#
# .. code-block:: cmake
#
#    find_package(Qt6 QUIET COMPONENTS OpenGLWidgets)
#
# and defines the variable :variable:`CGAL_Qt6_FOUND` and the function
# :command:`CGAL_setup_CGAL_Qt6_dependencies`.
#

if(CGAL_SetupCGAL_Qt6Dependencies_included)
  return()
endif()
set(CGAL_SetupCGAL_Qt6Dependencies_included TRUE)

#.rst:
# Used Modules
# ^^^^^^^^^^^^
#   - :module:`Qt6Config`

find_package(Qt6 QUIET COMPONENTS OpenGL OpenGLWidgets Widgets OPTIONAL_COMPONENTS Svg)

set(CGAL_Qt6_MISSING_DEPS "")
if(NOT Qt6OpenGLWidgets_FOUND)
  message( STATUS "NOTICE: NOT Qt6OpenGLWidgets_FOUND")
  set(CGAL_Qt6_MISSING_DEPS "Qt6OpenGLWidgets")
endif()
if(NOT Qt6_FOUND)
  message(STATUS "NOTICE: NOT Qt6_FOUND")
  set(CGAL_Qt6_MISSING_DEPS "${CGAL_Qt6_MISSING_DEPS} Qt6")
endif()
if(NOT EXISTS ${CGAL_GRAPHICSVIEW_PACKAGE_DIR}/include/CGAL/Qt/GraphicsItem.h)
  message(STATUS "NOTICE: NOT EXISTS GraphicsItem")
  set(CGAL_Qt6_MISSING_DEPS "${CGAL_Qt6_MISSING_DEPS} <CGAL/Qt/*.h> headers")
endif()

#.rst:
# Result Variables
# ^^^^^^^^^^^^^^^^
#
# .. variable:: CGAL_Qt6_FOUND
#
#    Set to `TRUE` if the dependencies of `CGAL_Qt6` were found.
#
if(NOT CGAL_Qt6_MISSING_DEPS)
  set(CGAL_Qt6_FOUND TRUE)
  set_property(GLOBAL PROPERTY CGAL_Qt6_FOUND TRUE)

  include(${CMAKE_CURRENT_LIST_DIR}/CGAL_Qt6_moc_and_resource_files.cmake)

  if(NOT TARGET CGAL_Qt6_moc_and_resources)
    add_library(CGAL_Qt6_moc_and_resources STATIC
      ${_CGAL_Qt6_MOC_FILES_private}
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
        ${_CGAL_Qt6_UI_FILES}
      ${_CGAL_Qt6_RESOURCE_FILES_private})
    target_include_directories( CGAL_Qt6_moc_and_resources PUBLIC ${CMAKE_CURRENT_BINARY_DIR})
    set_target_properties(CGAL_Qt6_moc_and_resources PROPERTIES
      POSITION_INDEPENDENT_CODE TRUE
      EXCLUDE_FROM_ALL TRUE
      AUTOMOC TRUE)
    target_link_libraries(CGAL_Qt6_moc_and_resources PUBLIC CGAL::CGAL Qt6::Widgets Qt6::OpenGLWidgets)
    if(Qt6Svg_FOUND)
      target_link_libraries(CGAL_Qt6_moc_and_resources PUBLIC Qt6::Svg)
    endif()
    add_library(CGAL::CGAL_Qt6_moc_and_resources ALIAS CGAL_Qt6_moc_and_resources)
    add_library(CGAL::Qt6_moc_and_resources ALIAS CGAL_Qt6_moc_and_resources)
  endif()

endif()

#get_property(QT_UIC_EXECUTABLE TARGET Qt6::uic PROPERTY LOCATION)
#message( STATUS "Qt6Core include:     ${Qt6Core_INCLUDE_DIRS}" )
#message( STATUS "Qt6 libraries:       ${Qt6Core_LIBRARIES} ${Qt6Gui_LIBRARIES} ${Qt6Svg_LIBRARIES} ${Qt6OpenGL_LIBRARIES}" )
#message( STATUS "Qt6Core definitions: ${Qt6Core_DEFINITIONS}" )
#message( STATUS "moc executable:      ${QT_MOC_EXECUTABLE}" )
#message( STATUS "uic executable:      ${QT_UIC_EXECUTABLE}" )

#.rst:
#
# Provided Functions
# ^^^^^^^^^^^^^^^^^^
#
# .. command:: CGAL_setup_CGAL_Qt6_dependencies
#
#   Link the target with the dependencies of `CGAL_Qt6`::
#
#     CGAL_setup_CGAL_Qt6_dependencies( target )
#
#   The dependencies are
#   added using :command:`target_link_libraries` with the ``INTERFACE``
#   keyword.
#
function(CGAL_setup_CGAL_Qt6_dependencies target)

  if($ENV{CGAL_FAKE_PUBLIC_RELEASE})
    target_compile_definitions( ${target} INTERFACE CGAL_FAKE_PUBLIC_RELEASE=1 )
  endif()
  target_link_libraries( ${target} INTERFACE CGAL::CGAL)
  target_link_libraries( ${target} INTERFACE CGAL::Qt6_moc_and_resources)
  target_link_libraries( ${target} INTERFACE Qt6::OpenGLWidgets )

  # Remove -Wdeprecated-copy, for g++ >= 9.0, because Qt6, as of
  # version 5.12, has a lot of [-Wdeprecated-copy] warnings.
  if( CMAKE_CXX_COMPILER_ID STREQUAL "GNU"
      AND NOT CMAKE_CXX_COMPILER_VERSION VERSION_LESS "9" )
    target_compile_options( ${target} INTERFACE "-Wno-deprecated-copy" "-Wno-cast-function-type" )
  endif()

endfunction()
