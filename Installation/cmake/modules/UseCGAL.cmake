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


  target_include_directories(CGAL INTERFACE "${CMAKE_CURRENT_BINARY_DIR}" )

  if(TARGET CGAL::CGAL)
    add_to_list( CGAL_LIBRARIES CGAL::CGAL )
  elseif(TARGET CGAL)
    add_to_list( CGAL_LIBRARIES CGAL )
  else()
    add_to_list( CGAL_LIBRARIES ${CGAL_LIBRARY} )
  endif()

  target_include_directories (CGAL INTERFACE ${CGAL_INCLUDE_DIRS})
  target_include_directories (CGAL SYSTEM INTERFACE ${CGAL_3RD_PARTY_INCLUDE_DIRS} )
  target_compile_definitions (CGAL INTERFACE ${CGAL_3RD_PARTY_DEFINITIONS}  ${CGAL_DEFINITIONS}  )
endif()
