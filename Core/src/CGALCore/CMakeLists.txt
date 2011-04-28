project (CGAL_Core)
   
# Minimal version of CMake:
cmake_minimum_required(VERSION 2.6.2)

# Tested version:
if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}" VERSION_GREATER 2.6)
  if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}.${CMAKE_PATCH_VERSION}" VERSION_GREATER 2.8.3)
    cmake_policy(VERSION 2.8.4)
  else()
    cmake_policy(VERSION 2.6)
  endif()
endif()

if ( NOT CGAL_FOUND )
  find_package(CGAL REQUIRED)
endif()

include(${CGAL_USE_FILE})

include_directories( ${CGAL_3RD_PARTY_INCLUDE_DIRS} )

build_cgal_library(CGAL_Core CGALCore "")

add_dependencies( CGAL_Core CGAL )

add_definitions ( ${CGAL_3RD_PARTY_DEFINITIONS} )
  
target_link_libraries( CGAL_Core ${CGAL_LIBRARY} ${CGAL_3RD_PARTY_LIBRARIES} )

set_target_properties(CGAL_Core PROPERTIES VERSION "${CGAL_SOVERSION}" SOVERSION "${CGAL_SONAME_VERSION}")

set( CGAL_Core_BASENAME CGAL_Core )

if (NOT CGAL_AUTO_LINK_ENABLED )

  set_target_properties( CGAL_Core PROPERTIES OUTPUT_NAME ${CGAL_Core_BASENAME} )
  target_link_libraries( CGAL_Core ${CGAL_3RD_PARTY_LIBRARIES})
  
  get_target_property( CGAL_Core_LIBRARY_FULL_PATH CGAL_Core LOCATION  )
  
  set(CGAL_Core_LIBRARY ${CGAL_Core_LIBRARY_FULL_PATH} CACHE STRING "The CGAL_Core library" FORCE)
  
  mark_as_advanced(CGAL_Core_LIBRARY)
  
else()

  set(CGAL_Core_LIBRARY)
  
  set(CGAL_Core_RELEASE_MANGLED_NAME        "${CGAL_Core_BASENAME}-${CGAL_TOOLSET}-mt" ) 
  set(CGAL_Core_DEBUG_MANGLED_NAME          "${CGAL_Core_BASENAME}-${CGAL_TOOLSET}-mt-gd" ) 
  set(CGAL_Core_MINSIZEREL_MANGLED_NAME     "${CGAL_Core_BASENAME}-${CGAL_TOOLSET}-mt-o" ) 
  set(CGAL_Core_RELWITHDEBINFO_MANGLED_NAME "${CGAL_Core_BASENAME}-${CGAL_TOOLSET}-mt-g" ) 
  
  set_target_properties( CGAL_Core PROPERTIES DEBUG_OUTPUT_NAME          "${CGAL_Core_DEBUG_MANGLED_NAME}" 
                                              RELEASE_OUTPUT_NAME        "${CGAL_Core_RELEASE_MANGLED_NAME}" 
                                              MINSIZEREL_OUTPUT_NAME     "${CGAL_Core_MINSIZEREL_MANGLED_NAME}" 
                                              RELWITHDEBINFO_OUTPUT_NAME "${CGAL_Core_RELWITHDEBINFO_MANGLED_NAME}" 
                      )
  
  if ( HAS_CFG_INTDIR )                      
    add_custom_command(TARGET CGAL_Core POST_BUILD COMMAND copy /Y \"$(TargetDir)CGAL_core*.lib\" \"$(SolutionDir)lib\" )
    add_custom_command(TARGET CGAL_Core POST_BUILD COMMAND if exist \"$(TargetDir)$(TargetName).pdb\" copy /Y \"$(TargetDir)CGAL_core*.pdb\" \"$(SolutionDir)lib\" )
  endif()
  
endif()

if ( CGAL_INSTALL_LIB_DIR )
  install(TARGETS CGAL_Core DESTINATION ${CGAL_INSTALL_LIB_DIR} )
endif()

