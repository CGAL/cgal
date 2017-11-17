if (NOT EXISTS ${CMAKE_BINARY_DIR}/Installation/include/CGAL/version.h)
  message(FATAL_ERROR "Cannot find Installation/include/CGAL/version.h. Make sure you are at the root of a CGAL branch")
endif()

file(READ "${CMAKE_BINARY_DIR}/Installation/include/CGAL/version.h" version_file_content)
string(REGEX MATCH "define CGAL_VERSION (.*)\n#define CGAL_VERSION_NR" CGAL_VERSION_FOUND "${version_file_content}")

if (CGAL_VERSION_FOUND)
  SET(CGAL_VERSION "${CMAKE_MATCH_1}")  
else()
  message(FATAL_ERROR "Cannot extract CGAL version number.")
endif()

if (NOT DEFINED DESTINATION)
  SET(DESTINATION "/tmp")
endif()

# Check consistency.
if(PUBLIC AND INTERNAL)
  message(FATAL_ERROR "You can only create either an internal or a public release.")
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
  if(IS_DIRECTORY ${pkg_dir} AND EXISTS ${pkg_dir}/package_info AND NOT "${pkg}" STREQUAL "Maintenance") # only consider packages
    if(VERBOSE)
      message(STATUS "handling ${pkg}")
    endif()

    # gather all files from this package
    set(all_files)
    file(GLOB_RECURSE pkg_files RELATIVE ${pkg_dir} ${pkg_dir}/*)
    # append the prefix
    foreach(f ${pkg_files})
      get_filename_component(fname ${f} NAME)
      if (NOT "${fname}" STREQUAL "TODO") # skip TODO files and unwanted directories
        list(APPEND all_files ${pkg_dir}/${f})
      endif()
    endforeach()

    # now copy them
    foreach(afile ${all_files})
      # get the path and remove the package part
      get_filename_component(afile_dir_tmp ${afile} PATH)
      string(REPLACE "${pkg_dir}" "" afile_dir ${afile_dir_tmp})
      if(NOT IS_DIRECTORY ${release_dir}/${afile_dir})
        file(MAKE_DIRECTORY ${release_dir}/${afile_dir})
      endif()

      file(COPY ${afile} DESTINATION ${release_dir}/${afile_dir})
    endforeach()
  endif()
endforeach()


file(REMOVE_RECURSE ${release_dir}/doc/fig_src)
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
