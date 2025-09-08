cmake_minimum_required(VERSION 3.16)
include_guard(GLOBAL)

function(process_CGAL_subdirectory entry subdir type_name)
  # For example, subdir can be "examples", type_name "example", and entry "Mesh_2"
  get_filename_component(ENTRY_DIR_NAME "${entry}" NAME)

  if( NOT "${CMAKE_SOURCE_DIR}" STREQUAL "${CMAKE_BINARY_DIR}") # out-of-source
    make_directory("${CMAKE_BINARY_DIR}/${subdir}/${ENTRY_DIR_NAME}")
  endif()

  message("-- Configuring ${subdir} in ${subdir}/${ENTRY_DIR_NAME}")
  list(APPEND CMAKE_MESSAGE_INDENT "  ")

  set(source_dir "")
  if(EXISTS ${entry}/CMakeLists.txt)
    set(source_dir ${entry})
  else()
    if(CGAL_CREATE_CMAKE_SCRIPT)
      execute_process(
        COMMAND bash ${CGAL_CREATE_CMAKE_SCRIPT} ${type_name} --source_dir "${entry}"
        WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/${subdir}/${ENTRY_DIR_NAME}"
        RESULT_VARIABLE RESULT_VAR OUTPUT_VARIABLE OUTPUT_VAR ERROR_VARIABLE ERROR_VAR)
      if(RESULT_VAR)
        message(AUTHOR_WARNING "Error with ${CGAL_CREATE_CMAKE_SCRIPT} ${type_name} --source_dir ${entry}\n${OUTPUT_VAR}\n${ERROR_VAR}")
      else()
        set(source_dir "${CMAKE_BINARY_DIR}/${subdir}/${ENTRY_DIR_NAME}")
      endif()
    endif()
  endif()
  if(source_dir)
    add_subdirectory( "${source_dir}" "${CMAKE_BINARY_DIR}/${subdir}/${ENTRY_DIR_NAME}" EXCLUDE_FROM_ALL)
  endif()
  list(POP_BACK CMAKE_MESSAGE_INDENT)
endfunction()

function(CGAL_handle_subdirectories subdir_name plural_name)
  if("${plural_name}" MATCHES "s$")
    string(LENGTH "${plural_name}" plural_name_length)
    math(EXPR plural_name_length_minus_one "${plural_name_length} - 1")
    string(SUBSTRING "${plural_name}" 0 "${plural_name_length_minus_one}" singular_name)
  else()
    set(singular_name "${plural_name}")
  endif()

  if(CGAL_BRANCH_BUILD)

    foreach(package ${CGAL_CONFIGURED_PACKAGES})
      #message (STATUS "Current package: ${package}")
      file(GLOB listtmp CONFIGURE_DEPENDS "${package}/${subdir_name}/*")
      list(APPEND list ${listtmp})
    endforeach()

  elseif(EXISTS "${CMAKE_SOURCE_DIR}/../../Installation/cmake/modules/CGAL_add_test.cmake")

    file(GLOB list CONFIGURE_DEPENDS "${CMAKE_SOURCE_DIR}/../../*/${subdir_name}/*")

  else()

    file(GLOB list CONFIGURE_DEPENDS "${subdir_name}/*")

  endif()

  if(NOT list)
    return()
  endif()

  list(SORT list)

  message("== Generating build files for ${plural_name} ==")
  foreach(entry ${list})

    if(IS_DIRECTORY ${entry})

      file(GLOB files "${entry}/*.cpp")

      # If there is no .cpp files, ignore the sub-directory
      if(files)
        process_CGAL_subdirectory("${entry}" ${subdir_name} ${singular_name})
      endif()

    endif()

  endforeach()
  message("== Generating build files for ${plural_name} (DONE) ==\n")

endfunction()
