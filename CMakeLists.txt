# Top level CMakeLists.txt for CGAL-branchbuild
project(CGAL CXX)

cmake_minimum_required(VERSION 2.6.2)

# option for branch build

option( CGAL_BRANCH_BUILD "Create CGAL from branch" ON)
mark_as_advanced( CGAL_BRANCH_BUILD )

# search for some SCM

# TODO

if ( EXISTS ${CMAKE_SOURCE_DIR}/.svn )
  set ( CGAL_SCM_NAME "svn" )
else ()
  if ( NOT EXISTS ${CMAKE_SOURCE_DIR}/.git )
     message( ERROR "Neither 'svn' nor 'git' as SCM found" ) 
  endif() 
  set ( CGAL_SCM_NAME "git" )
endif()

if ( ${CGAL_SCM_NAME} STREQUAL "svn" )

  find_program(SVN_EXECUTABLE svn DOC "subversion command line client")

  function(Subversion_GET_INFO dir variable)
    # use svnversion
    execute_process(COMMAND "${SVN_EXECUTABLE}" info --xml
      WORKING_DIRECTORY "${dir}"
      OUTPUT_VARIABLE ${variable}
      OUTPUT_STRIP_TRAILING_WHITESPACE
      RESULT_VARIABLE svn_error)
    if(svn_error)
      message("Warning: `svn info --xml \"${dir}\"` returned an error!")
      set(${variable} "")
    endif()
    set(${variable} ${${variable}} PARENT_SCOPE)
  endfunction()

  function(Subversion_GET_REVISION dir variable)
    Subversion_GET_INFO("${dir}" ${variable})
    string(REGEX REPLACE "^(.*\n)?   revision=\"([^\n]+)\".*url.*" "\\2" ${variable} "${${variable}}")
    if(NOT variable)
      message("Warning: can not get the Subversion revision of directory ${dir}")
    endif()
    set(${variable} ${${variable}} PARENT_SCOPE)
  endfunction(Subversion_GET_REVISION)

  function(Subversion_GET_URL dir variable)
    Subversion_GET_INFO("${dir}" ${variable})
    string(REGEX REPLACE ".*<url>([^<]+)</url>.*" "\\1" ${variable} "${${variable}}")
    if(NOT variable)
      message("Warning: can not get the Subversion URL of directory ${dir}")
    endif()
    set(${variable} ${${variable}} PARENT_SCOPE)
  endfunction(Subversion_GET_URL)

  Subversion_GET_URL("${CMAKE_CURRENT_SOURCE_DIR}" CGAL_TMP_SVN_BRANCH_NAME)
  # Remove the prefix of the URL "https://scm.gforge.inria.fr/".
  string(REGEX REPLACE "[a-z+]+://[a-z0-9.]+" "" CGAL_TMP_SVN_BRANCH_NAME "${CGAL_TMP_SVN_BRANCH_NAME}")

  set ( CGAL_SCM_BRANCH_NAME "${CGAL_TMP_SVN_BRANCH_NAME}" )
endif()

if ( ${CGAL_SCM_NAME} STREQUAL "git" )

  find_program(GIT_EXECUTABLE git DOC "git command line client")

#  EXECUTE_PROCESS(COMMAND ${GIT_EXECUTABLE} --git-dir=${CMAKE_SOURCE_DIR}/.git branch
#      OUTPUT_VARIABLE CGAL_GIT_BRANCH_OUT
#      OUTPUT_STRIP_TRAILING_WHITESPACE)
#  #foreach (line IN ${CGAL_GIT_BRANCH_OUT})
#  #  message (STATUS "Line: ${line}")
#    if ( STRING( REGEX MATCH "* " CGAL_GIT_BRANCH_LINE
#      ${CGAL_GIT_BRANCH_OUT} ) )
#       string ( REGEX REPLACE "# On branch " "" CGAL_GIT_BRANCH ${CGAL_GIT_BRANCH_LINE})
#    endif()
#  #endforeach()

  message( STATUS "Git branch ${CGAL_GIT_BRANCH}") 

  # TODO detect name
  set ( CGAL_SCM_BRANCH_NAME "n/a" )

endif()

# add option for duplicate file detection
option( CGAL_REPORT_DUPLICATE_FILES "Switch on to start (naive) detection of duplicate source- and headerfiles in packages" OFF)

# and finally start actual build
add_subdirectory( Installation )

