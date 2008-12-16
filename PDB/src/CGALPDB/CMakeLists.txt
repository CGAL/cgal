project (CGAL_PDB)
   
CMAKE_MINIMUM_REQUIRED(VERSION 2.4.5 FATAL_ERROR)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)
  
if ( COMMAND cmake_policy )
  cmake_policy( SET CMP0003 NEW )  
endif()

if ( NOT CGAL_FOUND )   
  find_package(CGAL REQUIRED)
  include(${CGAL_USE_FILE})
endif()

include_directories( ${CGAL_3RD_PARTY_INCLUDE_DIRS} )

add_library(CGAL_PDB "all_files.cpp")

set_target_properties(CGAL_PDB PROPERTIES VERSION "${CGAL_SOVERSION}" SOVERSION "${CGAL_SONAME_VERSION}")

target_link_libraries( CGAL_PDB ${CGAL_LIBRARY} ${CGAL_3RD_PARTY_LIBRARIES} )
        
add_definitions ( ${CGAL_3RD_PARTY_DEFINITIONS} )

add_dependencies( CGAL_PDB CGAL )

set( CGAL_PDB_BASENAME CGAL_PDB )

if (NOT CGAL_AUTO_LINK_ENABLED )

  set_target_properties( CGAL_PDB PROPERTIES OUTPUT_NAME "${CGAL_PDB_BASENAME}" )
  
  get_target_property( CGAL_PDB_LIBRARY_FULL_PATH CGAL_PDB LOCATION  )
  
  set(CGAL_PDB_LIBRARY "${CGAL_PDB_LIBRARY_FULL_PATH}" CACHE STRING "The CGAL_PDB library" FORCE)
  
  mark_as_advanced(CGAL_PDB_LIBRARY)
  
else()

  set(CGAL_PDB_LIBRARY)
  
  set(CGAL_PDB_RELEASE_MANGLED_NAME        "${CGAL_PDB_BASENAME}-${CGAL_TOOLSET}-mt" ) 
  set(CGAL_PDB_DEBUG_MANGLED_NAME          "${CGAL_PDB_BASENAME}-${CGAL_TOOLSET}-mt-gd" ) 
  set(CGAL_PDB_MINSIZEREL_MANGLED_NAME     "${CGAL_PDB_BASENAME}-${CGAL_TOOLSET}-mt-o" ) 
  set(CGAL_PDB_RELWITHDEBINFO_MANGLED_NAME "${CGAL_PDB_BASENAME}-${CGAL_TOOLSET}-mt-g" ) 
  
  set_target_properties( CGAL_PDB PROPERTIES DEBUG_OUTPUT_NAME          "${CGAL_PDB_DEBUG_MANGLED_NAME}" 
                                             RELEASE_OUTPUT_NAME        "${CGAL_PDB_RELEASE_MANGLED_NAME}" 
                                             MINSIZEREL_OUTPUT_NAME     "${CGAL_PDB_MINSIZEREL_MANGLED_NAME}" 
                                             RELWITHDEBINFO_OUTPUT_NAME "${CGAL_PDB_RELWITHDEBINFO_MANGLED_NAME}" 
                       )
  
  if ( HAS_CFG_INTDIR )                      
    add_custom_command(TARGET CGAL_PDB POST_BUILD COMMAND copy /Y \"$(TargetPath)\" \"$(SolutionDir)lib\" )
    add_custom_command(TARGET CGAL_PDB POST_BUILD COMMAND if exist \"$(TargetDir)$(TargetName).pdb\" copy /Y \"$(TargetDir)$(TargetName).pdb\" \"$(SolutionDir)lib\" )
  endif()
  
endif()

if ( CGAL_INSTALL_LIB_DIR )
  install(TARGETS CGAL_PDB DESTINATION ${CGAL_INSTALL_LIB_DIR} )
endif()
