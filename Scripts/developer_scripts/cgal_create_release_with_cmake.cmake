#option :
# DESTINATION the path where the release is created, default is /tmp
# PUBLIC=[ON/OFF] indicates if a public release should be built, default is OFF
# VERBOSE=[ON/OFF] makes the script more verbose, default is OFF
# CGAL_VERSION=release id used to update version.h, VERSION and the release directory. Can be 4.12-Ic-33, 4.12-I-32, 4.12, ... 
#   Must be followed by -beta<beta_number> if the release is a beta.
# CGAL_VERSION_NR=release string used to update version.h. Must be something like 1041200033 , or 10412009<beta number>0
# TESTSUITE=indicate if the release is meant to be used by the testsuite, default if OFF

if (NOT EXISTS ${CMAKE_BINARY_DIR}/Installation/include/CGAL/version.h)
  message(FATAL_ERROR "Cannot find Installation/include/CGAL/version.h. Make sure you are at the root of a CGAL branch")
endif()

file(READ "${CMAKE_BINARY_DIR}/Installation/include/CGAL/version.h" version_file_content)
string(REGEX MATCH "define CGAL_VERSION (.*)\n#define CGAL_VERSION_NR" CGAL_VERSION_FOUND "${version_file_content}")

if (CGAL_VERSION_FOUND)
  set(CGAL_VERSION_INPUT "${CMAKE_MATCH_1}")
  if (NOT CGAL_VERSION)
    set(CGAL_VERSION "${CGAL_VERSION_INPUT}")
  endif()
  set (GITHUB_PREFIX "https://github.com/CGAL/cgal/blob/releases/CGAL-${CGAL_VERSION}")
else()
  message(FATAL_ERROR "Cannot extract CGAL version number.")
endif()


if (NOT DEFINED DESTINATION)
  SET(DESTINATION "/tmp")
endif()

set(release_dir "${DESTINATION}/CGAL-${CGAL_VERSION}")
if(EXISTS ${release_dir})
  file(REMOVE_RECURSE ${release_dir})
endif()

if (PUBLIC)
  message(STATUS "Creating a public release ${CGAL_VERSION} in ${release_dir}")
else()
  message(STATUS "Creating an internal release ${CGAL_VERSION} in ${release_dir}")
endif()

file(MAKE_DIRECTORY "${release_dir}")
file(GLOB files RELATIVE ${CMAKE_BINARY_DIR} ${CMAKE_BINARY_DIR}/*)

foreach(pkg ${files})
  set(pkg_dir ${CMAKE_BINARY_DIR}/${pkg}) # use absolute path
  if(IS_DIRECTORY ${pkg_dir} AND (NOT "${pkg}" STREQUAL "Maintenance")
      AND (EXISTS ${pkg_dir}/package_info
           OR "${pkg}" STREQUAL "Documentation"
           OR "${pkg}" STREQUAL "Miscellany" ) ) # only consider packages
    if(VERBOSE)
      message(STATUS "handling ${pkg}")
    endif()

    # gather all files from this package
    set(all_files)
    file(GLOB_RECURSE pkg_files RELATIVE ${pkg_dir} ${pkg_dir}/*)
    # append the prefix
    foreach(f ${pkg_files})
      get_filename_component(fname ${f} NAME)
      if (NOT "${fname}" STREQUAL "TODO") # skip TODO files
        #make sure the target destination dir exists
        set(afile ${pkg_dir}/${f})
        get_filename_component(afile_dir_tmp ${afile} PATH)
        string(REPLACE "${pkg_dir}" "" afile_dir ${afile_dir_tmp})
        if(NOT IS_DIRECTORY ${release_dir}/${afile_dir})
          file(MAKE_DIRECTORY ${release_dir}/${afile_dir})
        endif()

        #copy the file (replace $URL$ and $ID$ for *.h and *.hpp)
        get_filename_component(fext ${fname} EXT)
        if ("${fext}" STREQUAL ".h" OR "${fext}" STREQUAL ".hpp")
          file(READ "${pkg_dir}/${f}" file_content)
          string(REPLACE "$URL$" "$URL: ${GITHUB_PREFIX}/${pkg}/${f} $" file_content "${file_content}")
          if(EXISTS ${CMAKE_BINARY_DIR}/.git)
            execute_process(
              COMMAND git --git-dir=${CMAKE_BINARY_DIR}/.git --work-tree=${CMAKE_BINARY_DIR} log -n1 "--format=format:%h %aI %an" -- "${pkg}/${f}"
              RESULT_VARIABLE RESULT_VAR
              OUTPUT_VARIABLE OUT_VAR
              )
            string(REPLACE "$Id$" "$Id: ${fname} ${OUT_VAR}" file_content "${file_content}")
          else()
            string(REPLACE "$Id$" "This file is from the release ${CGAL_VERSION} of CGAL" file_content "${file_content}")
          endif()
          file(WRITE ${release_dir}/${afile_dir}/${fname} "${file_content}")
        else()
          file(COPY ${afile} DESTINATION ${release_dir}/${afile_dir})
        endif()
      endif()
    endforeach()
    if (EXISTS "${release_dir}/doc/${pkg}")
      #generate filelist.txt used by doxygen ran on a release
      file(GLOB_RECURSE includes LIST_DIRECTORIES false RELATIVE "${CMAKE_BINARY_DIR}/${pkg}/include" "${CMAKE_BINARY_DIR}/${pkg}/include/CGAL/*.h")
      foreach(f ${includes})
        file(APPEND "${release_dir}/doc/${pkg}/filelist.txt" "${f}\n")
      endforeach()
      #remove fig_src directory
      if (IS_DIRECTORY "${release_dir}/doc/${pkg}/fig_src")
        file(REMOVE_RECURSE "${release_dir}/doc/${pkg}/fig_src")
      endif()
    endif()
  endif()
endforeach()

#create VERSION
file(WRITE ${release_dir}/VERSION "${CGAL_VERSION}")

#edit include/CGAL/version.h
file(READ "${release_dir}/include/CGAL/version.h" file_content)
#  update CGAL_GIT_HASH
if(EXISTS ${CMAKE_BINARY_DIR}/.git)
  execute_process(
    COMMAND git rev-parse HEAD
    RESULT_VARIABLE RESULT_VAR
    OUTPUT_VARIABLE OUT_VAR
    )
  string(REPLACE "CGAL_GIT_HASH abcdef\n" "CGAL_GIT_HASH ${OUT_VAR}" file_content "${file_content}")
endif()
#  update CGAL_RELEASE_DATE
string(TIMESTAMP TODAY "%Y%m%d")
string(REPLACE "CGAL_RELEASE_DATE 20170101" "CGAL_RELEASE_DATE ${TODAY}" file_content "${file_content}")
#  update CGAL_VERSION
string(REPLACE "CGAL_VERSION ${CGAL_VERSION_INPUT}" "CGAL_VERSION ${CGAL_VERSION}" file_content "${file_content}")
#  update CGAL_VERSION_NR
if (CGAL_VERSION_NR)
  string(REGEX REPLACE "CGAL_VERSION_NR 10[0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9]" "CGAL_VERSION_NR ${CGAL_VERSION_NR}" file_content "${file_content}")
endif()
file(WRITE ${release_dir}/include/CGAL/version.h "${file_content}")

# make an extra copy of examples and demos for the testsuite and generate
# create_cgal_test_with_cmake for tests, demos, and examples
if (TESTSUITE)
  file(GLOB tests RELATIVE "${release_dir}/test" "${release_dir}/test/*")
  foreach(d ${tests})
    if(IS_DIRECTORY "${release_dir}/test/${d}")
      if(NOT EXISTS "${release_dir}/test/${d}/cgal_test_with_cmake")
        execute_process(
          COMMAND ${CMAKE_BINARY_DIR}/Scripts/developer_scripts/create_cgal_test_with_cmake
          WORKING_DIRECTORY "${release_dir}/test/${d}"
          RESULT_VARIABLE RESULT_VAR
          OUTPUT_VARIABLE OUT_VAR
        )
        if(NOT "${RESULT_VAR}" STREQUAL "0")
          message(FATAL_ERROR "Error while running create_cgal_test_with_cmake in ${release_dir}/test/${d}")
        endif()
      endif()
    endif()
  endforeach()

  file(MAKE_DIRECTORY "${release_dir}/tmp")
  #copy demo/PKG to test/PKG_Demo
  file(GLOB demos RELATIVE "${release_dir}/demo" "${release_dir}/demo/*")
  foreach(d ${demos})
    if(IS_DIRECTORY "${release_dir}/demo/${d}")
      string(REGEX MATCH "^[a-z]" IS_RESOURCE_DIR ${d})
      if (IS_RESOURCE_DIR)
        file(COPY "${release_dir}/demo/${d}" DESTINATION "${release_dir}/test")
      else()
        #do the copy in 2 pass since we cannot specify the target name
        file(COPY "${release_dir}/demo/${d}" DESTINATION "${release_dir}/tmp")
        file(RENAME "${release_dir}/tmp/${d}" "${release_dir}/test/${d}_Demo")
        if(NOT EXISTS "${release_dir}/test/${d}_Demo/cgal_test_with_cmake")
          execute_process(
            COMMAND ${CMAKE_BINARY_DIR}/Scripts/developer_scripts/create_cgal_test_with_cmake --no-run
            WORKING_DIRECTORY "${release_dir}/test/${d}_Demo"
            RESULT_VARIABLE RESULT_VAR
            OUTPUT_VARIABLE OUT_VAR
          )
          if(NOT "${RESULT_VAR}" STREQUAL "0")
            message(FATAL_ERROR "Error while running create_cgal_test_with_cmake in ${release_dir}/test/${d}_Demo")
          endif()
        endif()
      endif()
    endif()
  endforeach()
  #copy examples/PKG to test/PKG_Examples
  file(GLOB examples RELATIVE "${release_dir}/examples" "${release_dir}/examples/*")
  foreach(d ${examples})
    if(IS_DIRECTORY "${release_dir}/examples/${d}")
      #do the copy in 2 pass since we cannot specify the target name
      file(COPY "${release_dir}/examples/${d}" DESTINATION "${release_dir}/tmp")
      file(RENAME "${release_dir}/tmp/${d}" "${release_dir}/test/${d}_Examples")
      if(NOT EXISTS "${release_dir}/test/${d}_Examples/cgal_test_with_cmake")
        execute_process(
          COMMAND ${CMAKE_BINARY_DIR}/Scripts/developer_scripts/create_cgal_test_with_cmake
          WORKING_DIRECTORY "${release_dir}/test/${d}_Examples"
          RESULT_VARIABLE RESULT_VAR
          OUTPUT_VARIABLE OUT_VAR
        )
        if(NOT "${RESULT_VAR}" STREQUAL "0")
          message(FATAL_ERROR "Error while running create_cgal_test_with_cmake in ${release_dir}/test/${d}_Examples")
        endif()
      endif()
    endif()
  endforeach()
  file(REMOVE_RECURSE "${release_dir}/tmp")
endif()

# removal of extra directories and files
file(REMOVE_RECURSE ${release_dir}/benchmark)
file(REMOVE_RECURSE ${release_dir}/archive)
file(REMOVE ${release_dir}/include/CGAL/license/generate_files.cmake)
file(REMOVE ${release_dir}/include/CGAL/license/README.md)
file(REMOVE ${release_dir}/include/CGAL/license/gpl.h.in)
file(REMOVE ${release_dir}/include/CGAL/license/package_list.txt)

if(PUBLIC) # we are not creating an internal release.
  # Taken from create_new_release.
  file(REMOVE_RECURSE ${release_dir}/test)
  file(REMOVE_RECURSE ${release_dir}/package_info)
  file(REMOVE_RECURSE ${release_dir}/developer_scripts)
  file(REMOVE_RECURSE ${release_dir}/doc)
  file(REMOVE_RECURSE ${release_dir}/include/CGAL/Test)
  file(REMOVE_RECURSE ${release_dir}/include/CGAL/Testsuite/)

  file(GLOB_RECURSE to_deletes RELATIVE ${release_dir} ${release_dir}/examples/*/cgal_test*)
  foreach(to_delete ${to_deletes})
    file(REMOVE ${release_dir}/${to_delete})
  endforeach()
  file(GLOB_RECURSE to_deletes RELATIVE ${release_dir} ${release_dir}/demo/*/cgal_test*)
  foreach(to_delete ${to_deletes})
    file(REMOVE ${release_dir}/${to_delete})
  endforeach()
endif()
