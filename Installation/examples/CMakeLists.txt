cmake_minimum_required(VERSION 3.1...3.15)
project(CGAL_EXAMPLES)

if(CGAL_BRANCH_BUILD)

  foreach(package ${CGAL_CONFIGURED_PACKAGES})
    #message (STATUS "Current package: ${package}")
    file(GLOB listtmp "${package}/examples/*")
    list(APPEND list ${listtmp})
  endforeach()

else()

  file(GLOB list "*")

endif()

list(SORT list)

if(NOT CGAL_BUILDING_LIBS)
  find_package(CGAL REQUIRED)
  include(${CGAL_MODULES_DIR}/CGAL_Macros.cmake)
endif()

message("== Generating build files for examples ==")
foreach(entry ${list})

  if(NOT ${entry} MATCHES ".*\\.svn\$" AND IS_DIRECTORY ${entry})

    file(GLOB files "${entry}/*.cpp")

    # If there is no .cpp files, ignore the sub-directory
    if(files)
      process_cgal_subdirectory("${entry}" examples example)
      # Note: process_CGAL_subdirectory is defined in cmake/modules/CGAL_Macros.cmake
    endif()

  endif()

endforeach()
message("== Generating build files for examples (DONE) ==\n")
