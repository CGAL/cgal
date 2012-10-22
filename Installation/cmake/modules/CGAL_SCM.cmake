# Common settings for CGAL cmake scripts
if( NOT CGAL_SCM_FILE_INCLUDED )
  set(CGAL_SCM_FILE_INCLUDED 1 )

if ( EXISTS ${CMAKE_SOURCE_DIR}/.git )
  set ( CGAL_SCM_NAME "git" )
else()
  if ( EXISTS ${CMAKE_SOURCE_DIR}/.svn )
    set ( CGAL_SCM_NAME "svn" )
  else ()
    message( ERROR "Neither 'svn' nor 'git' as SCM found" ) 
  endif()
endif()

if ( ${CGAL_SCM_NAME} STREQUAL "git" )

  find_program(GIT_EXECUTABLE git DOC "git command line client")
EXECUTE_PROCESS(COMMAND ${GIT_EXECUTABLE} --git-dir=${CMAKE_SOURCE_DIR}/.git symbolic-ref HEAD 
     OUTPUT_VARIABLE CGAL_GIT_BRANCH_OUT
     OUTPUT_STRIP_TRAILING_WHITESPACE)

  string ( REGEX REPLACE "refs/heads/" "" CGAL_GIT_BRANCH ${CGAL_GIT_BRANCH_OUT})
  message( STATUS "Git branch ${CGAL_GIT_BRANCH}") 

 set ( CGAL_SCM_BRANCH_NAME "${CGAL_GIT_BRANCH}" )

endif()
  
endif()
