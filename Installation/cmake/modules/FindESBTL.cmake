#This modules tries to find ESBTL
# Once done this will define
#
#  ESBTL_FOUND - system has ESBTL
#  ESBTL_INCLUDE_DIR - ESBTL include directory
#

# Is it already configured?
if (ESBTL_INCLUDE_DIR)
  set(ESBTL_FOUND TRUE)
else()  

  find_path(ESBTL_INCLUDE_DIR 
            NAMES ESBTL/default.h
            PATHS ENV ESBTL_INC_DIR
                  /usr/include
                  /usr/local/include
            DOC "The directory containing the ESBTL header files WITHOUT the ESBTL prefix"      
           )

  if(ESBTL_INCLUDE_DIR)
     set(ESBTL_FOUND TRUE)
  endif()
endif()

 
if(ESBTL_FOUND)
    message(STATUS "Found ESBTL: ${ESBTL_INCLUDE_DIR}")
    set(ESBTL_USE_FILE "UseESBTL")
endif()
