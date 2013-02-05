# Common settings for CGAL cmake scripts
if( NOT CGAL_SCM_FILE_INCLUDED )
  set(CGAL_SCM_FILE_INCLUDED 1 )

# CGAL_SCM_NAME can be either git or n/a as an indicator if we have a
# build system or not

# TODO: make that scm agnostic and turn it into a boolean
set(GIT_PARENT_DIR "${CMAKE_SOURCE_DIR}")
set(GIT_DIR "${GIT_PARENT_DIR}/.git")
if (EXISTS "${GIT_DIR}")
  set( CGAL_SCM_NAME "git" )
else()
  while(NOT EXISTS "${GIT_DIR}")	# .git dir not found, search parent directories
    set(GIT_PREVIOUS_PARENT "${GIT_PARENT_DIR}")
    get_filename_component(GIT_PARENT_DIR ${GIT_PARENT_DIR} PATH)
    if( "${GIT_PARENT_DIR}" STREQUAL "${GIT_PREVIOUS_PARENT}" )
      # We have reached the root directory, we are not in git
      set( _refspecvar "GITDIR-NOTFOUND" )
      set( _hashvar "GITDIR-NOTFOUND" )
      set( CGAL_SCM_NAME "n/a" )
      break()
    else()
      set( GIT_DIR "${GIT_PARENT_DIR}/.git" )
      set( CGAL_SCM_NAME "git" )
    endif()
  endwhile()
endif()

if ( "${CGAL_SCM_NAME}" STREQUAL "git" )
  find_program(GIT_EXECUTABLE git DOC "git command line client")
  execute_process(COMMAND ${GIT_EXECUTABLE} --git-dir=${GIT_DIR} rev-parse --symbolic --abbrev-ref HEAD
    OUTPUT_VARIABLE CGAL_GIT_BRANCH
    OUTPUT_STRIP_TRAILING_WHITESPACE)

  set( CGAL_SCM_BRANCH_NAME "${CGAL_GIT_BRANCH}" )

  execute_process(COMMAND ${GIT_EXECUTABLE} --git-dir=${GIT_DIR} rev-parse HEAD
    OUTPUT_VARIABLE CGAL_GIT_HASH
    OUTPUT_STRIP_TRAILING_WHITESPACE)

  # backward compatible
  set( CGAL_CREATED_SVN_REVISION "99999" )
else()
  set( CGAL_SCM_BRANCH_NAME "n/a")
  set( CGAL_CREATED_GIT_HASH "n/a" )
  set( CGAL_CREATED_SVN_REVISION "n/a" )
endif()

endif()
