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
  
  set( CGAL_LIBRARIES )

  foreach ( CGAL_COMPONENT ${CGAL_REQUESTED_COMPONENTS} )
    if(WITH_CGAL_${CGAL_COMPONENT})
      add_to_list( CGAL_LIBRARIES            ${CGAL_${CGAL_COMPONENT}_LIBRARY}              )
      add_to_list( CGAL_3RD_PARTY_LIBRARIES  ${CGAL_${CGAL_COMPONENT}_3RD_PARTY_LIBRARIES}  )
      
      add_to_list( CGAL_3RD_PARTY_INCLUDE_DIRS   ${CGAL_${CGAL_COMPONENT}_3RD_PARTY_INCLUDE_DIRS}   )
      add_to_list( CGAL_3RD_PARTY_DEFINITIONS    ${CGAL_${CGAL_COMPONENT}_3RD_PARTY_DEFINITIONS}    )
      add_to_list( CGAL_3RD_PARTY_LIBRARIES_DIRS ${CGAL_${CGAL_COMPONENT}_3RD_PARTY_LIBRARIES_DIRS} )
    endif()
  endforeach()
    
  include_directories( "${CMAKE_CURRENT_BINARY_DIR}" ) 

  add_to_list( CGAL_LIBRARIES ${CGAL_LIBRARY} )
  
  include_directories ( ${CGAL_3RD_PARTY_INCLUDE_DIRS} ${CGAL_INCLUDE_DIRS} )     
  add_definitions     ( ${CGAL_3RD_PARTY_DEFINITIONS}  ${CGAL_DEFINITIONS}  )
  
  link_directories    ( ${CGAL_LIBRARIES_DIR} ${CGAL_3RD_PARTY_LIBRARIES_DIRS} )
  link_libraries      ( ${CGAL_LIBRARIES}     ${CGAL_3RD_PARTY_LIBRARIES}      )

  
endif()
