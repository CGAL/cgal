macro (build_cgal_library LIBRARY_NAME LIBRARY_DIR_NAME ADDITIONAL_FILES)
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

# message(STATUS "Additional files: ${ADDITIONAL_FILES}")
  add_library (${LIBRARY_NAME} ${CMAKE_CURRENT_BINARY_DIR}/all_files.cpp ${ADDITIONAL_FILES})
#  add_library (${LIBRARY_NAME} ${CGAL_LIBRARY_SOURCE_FILES} ${ADDITIONAL_FILES}) # builing not creating temporary all_files.cpp

endmacro()

macro( configure_component DIR COMPONENT )
  option( WITH_CGAL_${COMPONENT} "Enable CGAL component ${COMPONENT}" ON)
  if ( WITH_CGAL_${COMPONENT})
    add_subdirectory( ${DIR} ${COMPONENT} )
  endif()
endmacro()

if(COMMAND unset)
  # Workaround: WITH_CGAL_Core can have been loaded from an old reference
  # cache Then remove it from the cache, and re-add it only if CGAL_Core
  # can be configured.
  unset(WITH_CGAL_Core CACHE)
endif()

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

      message(STATUS "Library config detected: ${CGAL_CONFIGURED_LIBRARY_NAME}")
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
   message(FATAL_ERROR "No library in configured packages found. Please fix package-file.")
endif()

list(REMOVE_DUPLICATES CGAL_CONFIGURED_LIBRARIES)
list(SORT CGAL_CONFIGURED_LIBRARIES)
cache_set(CGAL_CONFIGURED_LIBRARIES "${CGAL_CONFIGURED_LIBRARIES}")

message(STATUS "Library configs detected: ${CGAL_CONFIGURED_LIBRARIES}")

