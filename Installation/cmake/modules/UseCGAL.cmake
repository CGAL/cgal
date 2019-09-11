#
# UseCGAL.cmake can be included in a project to set the needed compiler and linker
# settings to use CGAL in a program.
#
# The variables used here are defined in the CGALConfig.cmake generated when CGAL was installed.
#
#
include(${CGAL_MODULES_DIR}/CGAL_Macros.cmake)

cgal_setup_module_path()

if(NOT USE_CGAL_FILE_INCLUDED)
  set(USE_CGAL_FILE_INCLUDED 1)

  include(${CMAKE_CURRENT_LIST_DIR}/CGAL_Common.cmake)
  if( CGAL_DEV_MODE OR RUNNING_CGAL_AUTO_TEST )
    include(${CMAKE_CURRENT_LIST_DIR}/CGAL_SetupFlags.cmake)
  else()
    include(${CMAKE_CURRENT_LIST_DIR}/CGAL_display_flags.cmake)
  endif()
  include(${CMAKE_CURRENT_LIST_DIR}/CGAL_GeneratorSpecificSettings.cmake)
  include(${CMAKE_CURRENT_LIST_DIR}/CGAL_TweakFindBoost.cmake)

  set( CGAL_LIBRARIES )

  foreach ( component ${CGAL_REQUESTED_COMPONENTS} )
    use_component( ${component} )
  endforeach()

  use_essential_libs()

  include_directories( "${CMAKE_CURRENT_BINARY_DIR}" )

  if(TARGET CGAL::CGAL)
    add_to_list( CGAL_LIBRARIES CGAL::CGAL )
  elseif(TARGET CGAL)
    add_to_list( CGAL_LIBRARIES CGAL )
  else()
    add_to_list( CGAL_LIBRARIES ${CGAL_LIBRARY} )
  endif()

  #message (STATUS "LIB: ${CGAL_LIBRARY}")
  #message (STATUS "LIBS: ${CGAL_LIBRARIES}")

  include_directories ( ${CGAL_INCLUDE_DIRS})
  include_directories ( SYSTEM ${CGAL_3RD_PARTY_INCLUDE_DIRS} )
  add_definitions     ( ${CGAL_3RD_PARTY_DEFINITIONS}  ${CGAL_DEFINITIONS}  )

  if (CGAL_HEADER_ONLY)
    if(NOT CGAL_NO_BLANKET_LINKING)
      link_directories    ( ${CGAL_3RD_PARTY_LIBRARIES_DIRS} )
      link_libraries      ( ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES}      )
    endif()
  else()
    if(NOT CGAL_NO_BLANKET_LINKING)
      link_directories    ( ${CGAL_LIBRARIES_DIR} ${CGAL_3RD_PARTY_LIBRARIES_DIRS} )
      link_libraries      ( ${CGAL_LIBRARIES}     ${CGAL_3RD_PARTY_LIBRARIES}      )
    endif()
  endif()

endif()
