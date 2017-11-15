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
#    find_package(OpenGL QUIET)
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
#   - :module:`FindOpenGL`
find_package(Qt5 QUIET COMPONENTS OpenGL Svg)
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
endif()

#get_property(QT_UIC_EXECUTABLE TARGET Qt5::uic PROPERTY LOCATION)
#message( STATUS "OpenGL include:      ${OPENGL_INCLUDE_DIR}" )
#message( STATUS "OpenGL libraries:    ${OPENGL_LIBRARIES}" )
#message( STATUS "OpenGL definitions:  ${OPENGL_DEFINITIONS}" )
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
  target_include_directories( ${target} SYSTEM ${keyword} ${OPENGL_INCLUDE_DIR})
  target_link_libraries( ${target} ${keyword} CGAL::CGAL)
  target_link_libraries( ${target} ${keyword} Qt5::OpenGL Qt5::Svg ${OPENGL_LIBRARIES})
endfunction()

include(${CMAKE_CURRENT_LIST_DIR}/Use_CGAL_Qt5_headers.cmake)
