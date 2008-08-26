#
# UseCGAL.cmake can be included in a project to set the needed compiler and linker
# settings to use CGAL in a program.
#
# The variables used here are defined in the CGALConfig.cmake generated when CGAL was installed.
#
#

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CGAL_CMAKE_MODULE_PATH})

if(NOT USE_CGAL_FILE_INCLUDED) 
  set(USE_CGAL_FILE_INCLUDED 1)

  include(CGALcommon)
  
  if ( IS_TOP_LEVEL )
    include(CGAL_SetupFlags)
    include(GeneratorSpecificSettings)
  endif()
  include_directories (${CGAL_BINARY_DIR}/include) # Plaform-specific include folder where compiler_config.h is located
  include_directories (${CGAL_INCLUDE_DIRS})             
  include_directories (${CGAL_3RD_PARTY_INCLUDE_DIRS})  

  add_definitions(${CGAL_3RD_PARTY_DEFINITIONS})
  
  link_directories( ${CGAL_LIBRARIES_DIR} ${CGAL_3RD_PARTY_LIBRARIES_DIRS} )
  
endif()
