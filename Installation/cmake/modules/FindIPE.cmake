# - Try to find Ipe
# Once done this will define
#
#  IPE_FOUND - system has Ipe
#  IPE_INCLUDE_DIR - the Ipe include directory
#  IPE_LIBRARIES - Link these to use Ipe
#


# Is it already configured?
if (IPE_INCLUDE_DIR AND IPE_LIBRARIES )
  set(IPE_FOUND TRUE)
else()  
  find_path(IPE_INCLUDE_DIR 
            NAMES ipelib.h
            PATHS /usr/include
                  /usr/local/include
           )

  find_library(IPE_LIBRARIES 
               NAMES ipe
               PATHS /usr/lib
                     /usr/local/lib
                     /usr/lib64
              )

  if(IPE_INCLUDE_DIR AND IPE_LIBRARIES)
    set(IPE_FOUND TRUE)
  endif()
endif()

if(IPE_FOUND)
    message(STATUS "Found Ipe: ${IPE_INCLUDE_DIR} ${IPE_LIBRARIES}")
endif()
