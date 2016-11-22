# - Try to find Ipe
# Once done this will define
#
#  IPE_FOUND - system has Ipe
#  IPE_INCLUDE_DIR - the Ipe include directory
#  IPE_LIBRARIES - Link these to use Ipe
#


# Is it already configured?
if (IPE_INCLUDE_DIR AND IPE_LIBRARIES AND IPE_FULL_VERSION)
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

  if(IPE_INCLUDE_DIR)
    file(READ "${IPE_INCLUDE_DIR}/ipebase.h" IPEBASE_H)
    string(REGEX MATCH "IPELIB_VERSION[ ]*=[ ]*([67])([0-9][0-9])([0-9][0-9]);" FOUND_IPE_VERSION "${IPEBASE_H}")
    if (FOUND_IPE_VERSION)
      set(IPE_VERSION ${CMAKE_MATCH_1} CACHE INTERNAL "Ipe version major number")
      set(IPE_MINOR_VERSION_1 ${CMAKE_MATCH_2} CACHE INTERNAL "Ipe version minor number")
      set(IPE_MINOR_VERSION_2 ${CMAKE_MATCH_3} CACHE INTERNAL "Ipe version patch number")
      set(IPE_FULL_VERSION "${IPE_VERSION}.${IPE_MINOR_VERSION_1}.${IPE_MINOR_VERSION_2}" CACHE INTERNAL "Ipe version x.y.z")
    endif()
  endif()

  if(IPE_INCLUDE_DIR AND IPE_LIBRARIES)
    set(IPE_FOUND TRUE)
  endif()
endif()

include(FindPackageHandleStandardArgs)

FIND_PACKAGE_HANDLE_STANDARD_ARGS(IPE
  REQUIRED_VARS IPE_INCLUDE_DIR IPE_LIBRARIES
  FOUND_VAR IPE_FOUND
  VERSION_VAR IPE_FULL_VERSION)

if(IPE_FOUND)
    message(STATUS "Found Ipe: ${IPE_INCLUDE_DIR} ${IPE_LIBRARIES}")
endif()
