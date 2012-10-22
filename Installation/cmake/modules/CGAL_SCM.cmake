# Common settings for CGAL cmake scripts
if( NOT CGAL_SCM_FILE_INCLUDED )
  set(CGAL_SCM_FILE_INCLUDED 1 )

set(GIT_PARENT_DIR "${CMAKE_SOURCE_DIR}")
set(GIT_DIR "${GIT_PARENT_DIR}/.git")
while(NOT EXISTS "${GIT_DIR}")	# .git dir not found, search parent directories
  set(GIT_PREVIOUS_PARENT "${GIT_PARENT_DIR}")
  get_filename_component(GIT_PARENT_DIR ${GIT_PARENT_DIR} PATH)
  if(GIT_PARENT_DIR STREQUAL GIT_PREVIOUS_PARENT)
    # We have reached the root directory, we are not in git
    set(${_refspecvar} "GITDIR-NOTFOUND" PARENT_SCOPE)
    set(${_hashvar} "GITDIR-NOTFOUND" PARENT_SCOPE)
  else()
    set(GIT_DIR "${GIT_PARENT_DIR}/.git")
  endif()
 endwhile()

if ( NOT "${_refspecvar}" STREQUAL "GITDIR-NOTFOUND" )
  set ( CGAL_SCM_NAME "git" )
endif()

#message ("GPD: ${GIT_PARENT_DIR}")

if ( ${CGAL_SCM_NAME} STREQUAL "git" )

  find_program(GIT_EXECUTABLE git DOC "git command line client")
  EXECUTE_PROCESS(COMMAND ${GIT_EXECUTABLE} --git-dir=${GIT_PARENT_DIR}/.git symbolic-ref HEAD 
     OUTPUT_VARIABLE CGAL_GIT_BRANCH_OUT
     OUTPUT_STRIP_TRAILING_WHITESPACE)

  string ( REGEX REPLACE "refs/heads/" "" CGAL_GIT_BRANCH ${CGAL_GIT_BRANCH_OUT})
  # message( STATUS "Git branch ${CGAL_GIT_BRANCH}") 

 set ( CGAL_SCM_BRANCH_NAME "${CGAL_GIT_BRANCH}" )

endif()
  
endif()
