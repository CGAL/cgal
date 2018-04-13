function(cgal_parse_version_h version_h_file name)
  file(READ "${version_h_file}" version_file_content)
  string(REGEX MATCH "define CGAL_VERSION (.*)\n#define *CGAL_VERSION_NR 1(0([0-9])|([1-9][0-9]))([0-9][0-9])([0-9])([0-9][0-9][0-9])" CGAL_VERSION_FOUND "${version_file_content}")
  if (NOT CGAL_VERSION_FOUND)
    message(FATAL_ERROR "Cannot extract CGAL version number.")
  endif()
  set(${name}  "${CMAKE_MATCH_1}" PARENT_SCOPE)

  # Now handle the optional arguments, to fill the version numbers
  # CMAKE_MATCH_2 corresponds to the alternative ([0-9])|([1-9][0-9])).
  # CMAKE_MATCH_3 and CMAKE_MATCH_4 corresponds to the two sub-expressions
  # of the alternative, and cannot be non-empty at the same time.
  set(${ARGV2} "${CMAKE_MATCH_3}${CMAKE_MATCH_4}" PARENT_SCOPE) # major version
  set(${ARGV3} "${CMAKE_MATCH_5}" PARENT_SCOPE)                 # minor version
  set(${ARGV4} "${CMAKE_MATCH_6}" PARENT_SCOPE)                 # patch number
  set(${ARGV5} "${CMAKE_MATCH_7}" PARENT_SCOPE)                 # build number
endfunction()
