function (collect_cgal_library LIBRARY_NAME ADDITIONAL_FILES)
  # IMPORTANT: First delete all_files.cpp
  if (EXISTS ${CMAKE_CURRENT_BINARY_DIR}/all_files.cpp)
    file(REMOVE ${CMAKE_CURRENT_BINARY_DIR}/all_files.cpp)
  endif()

  # THEN collect *.cpp
  foreach (package ${CGAL_CONFIGURED_PACKAGES} )
    file(GLOB CGAL_LIBRARY_SOURCE_FILES_TMP ${package}/src/${LIBRARY_NAME}/*.cpp)
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

  if(CGAL_AUTO_LINK_ENABLED)
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
  endif()

  install(TARGETS ${LIBRARY_NAME} EXPORT CGALExports
          RUNTIME DESTINATION ${CGAL_INSTALL_BIN_DIR}
          LIBRARY DESTINATION ${CGAL_INSTALL_LIB_DIR}
          ARCHIVE DESTINATION ${CGAL_INSTALL_LIB_DIR})
endfunction()

function( configure_component DIR COMPONENT )
  if(${COMPONENT} STREQUAL "CGAL_Qt3")
    set(enabled_by_default OFF)
  else()
    set(enabled_by_default ON)
  endif()
  option( WITH_${COMPONENT} "Enable CGAL component ${COMPONENT}" ${enabled_by_default})
  if ( WITH_${COMPONENT})
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

# TODO: Seems useless, because it is called again in ../CMakeLists.txt
#       Should probably be removed. -- Laurent Rineau, 2014/07/22
create_CGALconfig_files()

set(CGAL_DIR ${CMAKE_BINARY_DIR})

add_subdirectory(CGAL)

# search libs
set(CGAL_CONFIGURED_LIBRARIES "")
 
foreach(package ${CGAL_CONFIGURED_PACKAGES}) 
  file(GLOB CONFIGURED_LIBS_IN_PACKAGE ${package}/src/CGAL_*/CMakeLists.txt)
  foreach (libconfigfile ${CONFIGURED_LIBS_IN_PACKAGE})
    get_filename_component(CGAL_CONFIGURED_LIBRARY_NAME ${libconfigfile} PATH)
    get_filename_component(CGAL_CONFIGURED_LIBRARY_NAME ${CGAL_CONFIGURED_LIBRARY_NAME} NAME)

    if (NOT ${CGAL_CONFIGURED_LIBRARY_NAME} STREQUAL "CGAL")

     message(STATUS "Sources for CGAL component library '${CGAL_CONFIGURED_LIBRARY_NAME}' detected")
      list(APPEND CGAL_CONFIGURED_LIBRARIES ${CGAL_CONFIGURED_LIBRARY_NAME})
      #message(STATUS "Library config detected: ${CGAL_CONFIGURED_LIBRARIES}")

      cache_set(${CGAL_CONFIGURED_LIBRARY_NAME}_LIBRARY                  "" )
      cache_set(${CGAL_CONFIGURED_LIBRARY_NAME}_3RD_PARTY_DEFINITIONS    "" )
      cache_set(${CGAL_CONFIGURED_LIBRARY_NAME}_3RD_PARTY_INCLUDE_DIRS   "" )
      cache_set(${CGAL_CONFIGURED_LIBRARY_NAME}_3RD_PARTY_LIBRARIES      "" )
      cache_set(${CGAL_CONFIGURED_LIBRARY_NAME}_3RD_PARTY_LIBRARIES_DIRS "" ) 
      if (${CGAL_CONFIGURED_LIBRARY_NAME} STREQUAL "CGAL_Core") 
        if (NOT CGAL_NO_CORE)
          configure_component( ${package}/src/${CGAL_CONFIGURED_LIBRARY_NAME} ${CGAL_CONFIGURED_LIBRARY_NAME})
        else(NOT CGAL_NO_CORE)
          message( STATUS "CGAL_Core is not being configured (missing GMP, or 64bit architecture).")
        endif(NOT CGAL_NO_CORE)
      else()
        configure_component( ${package}/src/${CGAL_CONFIGURED_LIBRARY_NAME} ${CGAL_CONFIGURED_LIBRARY_NAME})
      endif()

     cache_get(${CGAL_CONFIGURED_LIBRARY_NAME}_LIBRARY)

     cache_get(${CGAL_CONFIGURED_LIBRARY_NAME}_3RD_PARTY_DEFINITIONS   )
     cache_get(${CGAL_CONFIGURED_LIBRARY_NAME}_3RD_PARTY_INCLUDE_DIRS  )
     cache_get(${CGAL_CONFIGURED_LIBRARY_NAME}_3RD_PARTY_LIBRARIES     )
     cache_get(${CGAL_CONFIGURED_LIBRARY_NAME}_3RD_PARTY_LIBRARIES_DIRS) 

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

# CGAL_CONFIGURED_LIBRARIES is a virtually useless variable, because
# we do not know if those components actually created a target.  It
# might seem bad that we hardcode the target names here, but they have
# always been hardcoded in the CGALConfig files, so we do not make it
# any worse. We can also not use export(EXPORT) since that is only
# available in newer CMake versions.
set(CGAL_ACTUAL_CONFIGURED_LIBRARIES "")
macro(CGAL_add_if_target TARGET)
  if(TARGET ${TARGET})
    list(APPEND CGAL_ACTUAL_CONFIGURED_LIBRARIES ${TARGET})
  endif()
endmacro()

CGAL_add_if_target(CGAL)
CGAL_add_if_target(CGAL_Core)
CGAL_add_if_target(CGAL_ImageIO)
CGAL_add_if_target(CGAL_Qt3)
CGAL_add_if_target(CGAL_Qt5)

export(TARGETS ${CGAL_ACTUAL_CONFIGURED_LIBRARIES} NAMESPACE CGAL:: FILE "${CMAKE_BINARY_DIR}/CGALExports.cmake")

message(STATUS "Sources for CGAL component libraries '${CGAL_CONFIGURED_LIBRARIES}' detected")

