function (collect_cgal_library LIBRARY_NAME LIBRARY_DIR_NAME ADDITIONAL_FILES)
  # IMPORTANT: First delete all_files.cpp
  if (EXISTS ${CMAKE_CURRENT_BINARY_DIR}/all_files.cpp)
    file(REMOVE ${CMAKE_CURRENT_BINARY_DIR}/all_files.cpp)
  endif()

  # THEN collect *.cpp
  foreach (package ${CGAL_CONFIGURED_PACKAGES} )
    file(GLOB CGAL_LIBRARY_SOURCE_FILES_TMP ${package}/src/${LIBRARY_DIR_NAME}/*.cpp)
    list(APPEND CGAL_LIBRARY_SOURCE_FILES ${CGAL_LIBRARY_SOURCE_FILES_TMP})
  endforeach()

  foreach(source ${CGAL_LIBRARY_SOURCE_FILES})
    # It may happen that a file all_files.cpp had been created in-source by
    # a previous in-source build. We certainly do not want to include that
    # file in the new all_files.cpp because .cpp files would be included
    # twice, and that breaks the one-definition rule (ODR).
    if(NOT source MATCHES ".*/all_files.cpp")
      file(APPEND ${CMAKE_CURRENT_BINARY_DIR}/all_files.cpp "\#include \"${source}\"\n")
    endif()
  endforeach()

  if(CGAL_BUILD_SHARED_LIBS)
    set(CGAL_LIB_PREFIX "")
  else()
    set(CGAL_LIB_PREFIX "lib")
  endif()

  set(RELEASE_MANGLED_NAME        "${CGAL_LIB_PREFIX}${LIBRARY_NAME}-${CGAL_TOOLSET}-mt-${CGAL_FULL_VERSION}" ) 
  set(DEBUG_MANGLED_NAME          "${CGAL_LIB_PREFIX}${LIBRARY_NAME}-${CGAL_TOOLSET}-mt-gd-${CGAL_FULL_VERSION}" ) 

  if(CGAL_BUILD_SHARED_LIBS AND WIN32)
    configure_file(${CGAL_INSTALLATION_PACKAGE_DIR}/src/CGAL_libs_verinfo.rc.in ${LIBRARY_NAME}_verinfo.rc @ONLY)
    set(rc_file ${CMAKE_CURRENT_BINARY_DIR}/${LIBRARY_NAME}_verinfo.rc)
  else()
    set(rc_file )
  endif()

  add_library (${LIBRARY_NAME}
    ${CMAKE_CURRENT_BINARY_DIR}/all_files.cpp 
    ${rc_file}
    ${ADDITIONAL_FILES})
#  add_library (${LIBRARY_NAME} ${CGAL_LIBRARY_SOURCE_FILES} ${rc_file} ${ADDITIONAL_FILES}) # builing not creating temporary all_files.cpp
  if(CGAL_SOVERSION AND CGAL_SONAME_VERSION)
    set_target_properties(${LIBRARY_NAME} PROPERTIES 
      VERSION "${CGAL_SOVERSION}" 
      SOVERSION "${CGAL_SONAME_VERSION}")
  endif()

  if (NOT CGAL_AUTO_LINK_ENABLED )
    get_target_property( ${LIBRARY_NAME}_LIBRARY_FULL_PATH ${LIBRARY_NAME} LOCATION  )
  
    set(${LIBRARY_NAME}_LIBRARY ${${LIBRARY_NAME}_LIBRARY_FULL_PATH} CACHE STRING "The ${LIBRARY_NAME} library" FORCE)
    mark_as_advanced(${LIBRARY_NAME}_LIBRARY)
  else(NOT CGAL_AUTO_LINK_ENABLED)
    set(${LIBRARY_NAME}_LIBRARY)

    set_target_properties( ${LIBRARY_NAME} PROPERTIES 
      OUTPUT_NAME_DEBUG          "${DEBUG_MANGLED_NAME}" 
      OUTPUT_NAME_RELEASE        "${RELEASE_MANGLED_NAME}" 
      OUTPUT_NAME_MINSIZEREL     "${RELEASE_MANGLED_NAME}" 
      OUTPUT_NAME_RELWITHDEBINFO "${RELEASE_MANGLED_NAME}" 
      )
    
    if ( HAS_CFG_INTDIR )
      add_custom_command(TARGET ${LIBRARY_NAME} POST_BUILD COMMAND if exist \"$(TargetDir)$(TargetName).dll\" copy /Y \"$(TargetDir)$(TargetName).dll\" \"$(TargetDir)..\" )
      add_custom_command(TARGET ${LIBRARY_NAME} POST_BUILD COMMAND if exist \"$(TargetDir)$(TargetName).pdb\" copy /Y \"$(TargetDir)$(TargetName).pdb\" \"$(TargetDir)..\" )
      add_custom_command(TARGET ${LIBRARY_NAME} POST_BUILD COMMAND if exist \"$(SolutionDir)lib\\$(ConfigurationName)\\$(TargetName).lib\" copy /Y \"$(SolutionDir)lib\\$(ConfigurationName)\\$(TargetName).lib\" \"$(SolutionDir)lib\" )
    endif()
  endif(NOT CGAL_AUTO_LINK_ENABLED)

  install(TARGETS ${LIBRARY_NAME}
          RUNTIME DESTINATION ${CGAL_INSTALL_BIN_DIR}
          LIBRARY DESTINATION ${CGAL_INSTALL_LIB_DIR}
          ARCHIVE DESTINATION ${CGAL_INSTALL_LIB_DIR})
endfunction()

function( configure_component DIR COMPONENT )
  option( WITH_CGAL_${COMPONENT} "Enable CGAL component ${COMPONENT}" ON)
  if ( WITH_CGAL_${COMPONENT})
    add_subdirectory( ${DIR} ${COMPONENT} )
  endif()
endfunction()

# Output directory of libraries (.so, .dylib, non-static .lib)
set(CGAL_LIBRARIES_DIR ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CGAL_LIBRARIES_DIR})
# Output directory of static libraries (.a, static .lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CGAL_LIBRARIES_DIR})

# For output directory of DLLs (.dll files)
set(CGAL_RUNTIME_DIR ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CGAL_RUNTIME_DIR})

## Prevent use of intermediate directory with MSVC, because we use name
## mangling
#foreach(Conf_type ${CMAKE_CONFIGURATION_TYPES})
#  # Conf_type is Debug, Release, DebWithDebInfo...
#  string(TOUPPER ${Conf_type} CONF_TYPE)
#  # CONF_TYPE is DEBUG, RELEASE, DEBWITHDEBINFO
#  set(CMAKE_LIBRARY_OUTPUT_DIRECTORY_${CONF_TYPE} ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
#  set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY_${CONF_TYPE} ${CMAKE_ARCHIVE_OUTPUT_DIRECTORY})
#  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_${CONF_TYPE} ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
#endforeach()

create_CGALconfig_files()

set(CGAL_DIR ${CMAKE_BINARY_DIR})

add_subdirectory(CGAL)

# search libs
set(CGAL_CONFIGURED_LIBRARIES "" PARENT_SCOPE)
 
foreach(package ${CGAL_CONFIGURED_PACKAGES}) 
  file(GLOB CONFIGURED_LIBS_IN_PACKAGE ${package}/src/CGAL*/CMakeLists.txt)
  foreach (libconfigfile ${CONFIGURED_LIBS_IN_PACKAGE})
    string(REPLACE "${package}/src/" "" libconfigfile ${libconfigfile})
    string(REPLACE "//CMakeLists.txt" "" CGAL_CONFIGURED_LIBRARY_SUBDIR ${libconfigfile})
    if (NOT ${CGAL_CONFIGURED_LIBRARY_SUBDIR} STREQUAL "CGAL")

      string(REPLACE "CGAL" "" CGAL_CONFIGURED_LIBRARY_NAME ${CGAL_CONFIGURED_LIBRARY_SUBDIR})

      message(STATUS "Sources for CGAL component library '${CGAL_CONFIGURED_LIBRARY_NAME}' detected")
      list(APPEND CGAL_CONFIGURED_LIBRARIES ${CGAL_CONFIGURED_LIBRARY_NAME})
      #message(STATUS "Library config detected: ${CGAL_CONFIGURED_LIBRARIES}")

      cache_set(CGAL_${CGAL_CONFIGURED_LIBRARY_NAME}_LIBRARY                  "" )

      cache_set(CGAL_${CGAL_CONFIGURED_LIBRARY_NAME}_3RD_PARTY_DEFINITIONS    "" )
      cache_set(CGAL_${CGAL_CONFIGURED_LIBRARY_NAME}_3RD_PARTY_INCLUDE_DIRS   "" )
      cache_set(CGAL_${CGAL_CONFIGURED_LIBRARY_NAME}_3RD_PARTY_LIBRARIES      "" )
      cache_set(CGAL_${CGAL_CONFIGURED_LIBRARY_NAME}_3RD_PARTY_LIBRARIES_DIRS "" ) 

      if (${CGAL_CONFIGURED_LIBRARY_NAME} STREQUAL "Core") 
        if (NOT CGAL_NO_CORE)
          configure_component( ${package}/src/${CGAL_CONFIGURED_LIBRARY_SUBDIR} ${CGAL_CONFIGURED_LIBRARY_NAME})
        else(NOT CGAL_NO_CORE)
          message( STATUS "CGAL_Core is not being configured (missing GMP, or 64bit architecture).")
        endif(NOT CGAL_NO_CORE)
      else()
        configure_component( ${package}/src/${CGAL_CONFIGURED_LIBRARY_SUBDIR} ${CGAL_CONFIGURED_LIBRARY_NAME})
      endif()

     cache_get(CGAL_${CGAL_CONFIGURED_LIBRARY_NAME}_LIBRARY)

     cache_get(CGAL_${CGAL_CONFIGURED_LIBRARY_NAME}_3RD_PARTY_DEFINITIONS   )
     cache_get(CGAL_${CGAL_CONFIGURED_LIBRARY_NAME}_3RD_PARTY_INCLUDE_DIRS  )
     cache_get(CGAL_${CGAL_CONFIGURED_LIBRARY_NAME}_3RD_PARTY_LIBRARIES     )
     cache_get(CGAL_${CGAL_CONFIGURED_LIBRARY_NAME}_3RD_PARTY_LIBRARIES_DIRS) 

    else() 
 
      list(APPEND CGAL_CONFIGURED_LIBRARIES CGAL)

    endif()
  endforeach()
endforeach()

if (NOT CGAL_CONFIGURED_LIBRARIES) 
   message(FATAL_ERROR "No component library in configured packages found. Please fix package-file.")
endif()

list(REMOVE_DUPLICATES CGAL_CONFIGURED_LIBRARIES)
list(SORT CGAL_CONFIGURED_LIBRARIES)
cache_set(CGAL_CONFIGURED_LIBRARIES "${CGAL_CONFIGURED_LIBRARIES}")

message(STATUS "Sources for CGAL component libraries '${CGAL_CONFIGURED_LIBRARIES}' detected")

