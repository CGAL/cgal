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
    message (STATUS "CGAL requested component: ${CGAL_COMPONENT}")

    if(WITH_CGAL_${CGAL_COMPONENT})
      add_to_list( CGAL_LIBRARIES            ${CGAL_${CGAL_COMPONENT}_LIBRARY}              )
      add_to_list( CGAL_3RD_PARTY_LIBRARIES  ${CGAL_${CGAL_COMPONENT}_3RD_PARTY_LIBRARIES}  )
      
      add_to_list( CGAL_3RD_PARTY_INCLUDE_DIRS   ${CGAL_${CGAL_COMPONENT}_3RD_PARTY_INCLUDE_DIRS}   )
      add_to_list( CGAL_3RD_PARTY_DEFINITIONS    ${CGAL_${CGAL_COMPONENT}_3RD_PARTY_DEFINITIONS}    )
      add_to_list( CGAL_3RD_PARTY_LIBRARIES_DIRS ${CGAL_${CGAL_COMPONENT}_3RD_PARTY_LIBRARIES_DIRS} )
    endif()

    # TODO-EBEB enable ALL_PRECONFIGURED_LIBS with cmake-option (which is, e.g., not given by testsuite
    if ( ${CGAL_COMPONENT} STREQUAL "ALL_PRECONFIGURED_LIBS" )

      message( STATUS "External libraries are all used")
      foreach ( CGAL_3RD_PARTY_LIB ${CGAL_SUPPORTING_3RD_PARTY_LIRARIES})
        if (${CGAL_3RD_PARTY_LIB}_FOUND) 
          use_lib( ${CGAL_3RD_PARTY_LIB} "###${${CGAL_3RD_PARTY_LIB}_USE_FILE}")
        endif()
      endforeach()
    
    else() 

      if ( ${CGAL_COMPONENT}_FOUND) 

        message( STATUS "External library ${CGAL_COMPONENT} is preconfigured")
        use_lib( ${CGAL_COMPONENT} "###${${CGAL_COMPONENT}_USE_FILE}")

      else()

        message( STATUS "External library ${CGAL_COMPONENT} is not preconfigured")
        find_package( ${CGAL_COMPONENT} )
        if (${CGAL_COMPONENT}_FOUND) 
          use_lib( ${CGAL_COMPONENT} "###${${CGAL_COMPONENT}_USE_FILE}")
        endif()
     
      endif()
    endif()
  endforeach()
    
  include_directories( "${CMAKE_CURRENT_BINARY_DIR}" ) 

  # need to get variable from cache while compiling CGAL, while in a demo it is set in CGALConfig.cmake
  if ( NOT CGAL_LIBRARY ) 
    cache_get(CGAL_LIBRARY)
  endif()
  add_to_list( CGAL_LIBRARIES ${CGAL_LIBRARY} )

  #message (STATUS "LIB: ${CGAL_LIBRARY}")
  #message (STATUS "LIBS: ${CGAL_LIBRARIES}")
  
  include_directories ( ${CGAL_3RD_PARTY_INCLUDE_DIRS} ${CGAL_INCLUDE_DIRS} )     
  add_definitions     ( ${CGAL_3RD_PARTY_DEFINITIONS}  ${CGAL_DEFINITIONS}  )
  
  link_directories    ( ${CGAL_LIBRARIES_DIR} ${CGAL_3RD_PARTY_LIBRARIES_DIRS} )
  link_libraries      ( ${CGAL_LIBRARIES}     ${CGAL_3RD_PARTY_LIBRARIES}      )

  
endif()
