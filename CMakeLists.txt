# Top level CMakeLists.txt for CGAL-branchbuild
CMAKE_MINIMUM_REQUIRED(VERSION 2.6.2)

# option for branch build

option( CGAL_BRANCH_BUILD "Create CGAL from branch" ON)

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

  FIND_PROGRAM(SVN_EXECUTABLE svn DOC "subversion command line client")

  # TODO detect name
  set ( CGAL_SCM_BRANCH_NAME "n/a" )

endif()

if ( ${CGAL_SCM_NAME} STREQUAL "git" )

  FIND_PROGRAM(GIT_EXECUTABLE git DOC "git command line client")

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

