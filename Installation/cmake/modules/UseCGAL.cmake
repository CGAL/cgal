#
# UseCGAL.cmake can be included in a project to set the needed compiler and linker
# settings to use CGAL in a program.
#
# The variables used here are defined in the CGALConfig.cmake generated when CGAL was installed.
#
#
set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

include(${CGAL_MODULES_DIR}/CGAL_Macros.cmake)

cgal_setup_module_path()

if(NOT USE_CGAL_FILE_INCLUDED) 
  set(USE_CGAL_FILE_INCLUDED 1)

  include(CGAL_Common)
  include(CGAL_SetupFlags)
  include(CGAL_GeneratorSpecificSettings)
  include(CGAL_TweakFindBoost)
  
  set( CGAL_LIBRARIES )

  foreach ( component ${CGAL_REQUESTED_COMPONENTS} )
    use_component( ${component} )
  endforeach()

  use_essential_libs()
    
  include_directories( "${CMAKE_CURRENT_BINARY_DIR}" ) 

  # need to get variable from cache while compiling CGAL, while in a demo it is set in CGALConfig.cmake
  if ( NOT CGAL_LIBRARY ) 
    cache_get(CGAL_LIBRARY)
  endif()
  add_to_list( CGAL_LIBRARIES ${CGAL_LIBRARY} )

  #message (STATUS "LIB: ${CGAL_LIBRARY}")
  #message (STATUS "LIBS: ${CGAL_LIBRARIES}")
  
  include_directories ( ${CGAL_INCLUDE_DIRS})
  include_directories ( SYSTEM ${CGAL_3RD_PARTY_INCLUDE_DIRS} )
  add_definitions     ( ${CGAL_3RD_PARTY_DEFINITIONS}  ${CGAL_DEFINITIONS}  )
  
  link_directories    ( ${CGAL_LIBRARIES_DIR} ${CGAL_3RD_PARTY_LIBRARIES_DIRS} )
  link_libraries      ( ${CGAL_LIBRARIES}     ${CGAL_3RD_PARTY_LIBRARIES}      )

  
endif()
