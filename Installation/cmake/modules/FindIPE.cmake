# - Try to find Ipe
# Once done this will define
#
#  IPE_FOUND - system has Ipe
#  IPE_INCLUDE_DIR - the Ipe include directory
#  IPE_LIBRARIES - Link these to use Ipe
#  WITH_IPE_7 - indicates if the compatibility with the version 7 of IPE must be used
#


macro( remove_leading_zero var )
  string(SUBSTRING "${${var}}" 0 1 ONECHAR)
  string(COMPARE EQUAL "${ONECHAR}" "0" ISZERO)
  if (${ISZERO})
    string(SUBSTRING "${${var}}" 1 1 ONECHAR)
    set(${var} ${ONECHAR})
  endif()
endmacro()

# Is it already configured?
if (IPE_INCLUDE_DIR AND IPE_LIBRARY_DIR )
  set(IPE_FOUND TRUE)
else()  

  find_path(IPE_INCLUDE_DIR 
            NAMES ipelib.h
            PATHS /usr/include
                  /usr/local/include
           )

  find_library(IPE_LIBRARY 
               NAMES ipe
               PATHS /usr/lib
                     /usr/local/lib
                     /usr/lib64
              )

  if(IPE_INCLUDE_DIR AND IPE_LIBRARY)
    set(IPE_FOUND TRUE)
    if (${IPE_VERSION} STREQUAL "AUTODETECT")
      FILE(READ "${IPE_INCLUDE_DIR}/ipebase.h" IPEBASE_H)
      STRING(REGEX MATCH "IPELIB_VERSION[ ]*=[ ]*([67])([0-9][0-9])([0-9][0-9]);" found_ipe_version "${IPEBASE_H}")
      if (found_ipe_version)
        set(IPE_VERSION ${CMAKE_MATCH_1})
        set(IPE_MINOR_VERSION_1 ${CMAKE_MATCH_2})
        set(IPE_MINOR_VERSION_2 ${CMAKE_MATCH_3})
      endif()
    endif()
    if (${IPE_VERSION} EQUAL "7")
      set(WITH_IPE_7 ON)
    elseif(${IPE_VERSION} EQUAL "6")
      set(WITH_IPE_7 OFF)
    else()
      message("-- Error: ${IPE_VERSION} is not a supported version of IPE (only 6 and 7 are).")
      set(IPE_FOUND FALSE)
    endif()        
  endif()
  get_filename_component(IPE_LIBRARY_DIR ${IPE_LIBRARY} PATH)
endif()

if (IPE_FOUND AND NOT IPELET_INSTALL_DIR)
  message("-- Using IPE version ${IPE_VERSION} compatibility.") 
  if (WITH_IPE_7)
    remove_leading_zero(IPE_MINOR_VERSION_1)
    remove_leading_zero(IPE_MINOR_VERSION_2)
    set(INSTALL_PATHS ${INSTALL_PATHS} "${IPE_LIBRARY_DIR}/ipe/7.${IPE_MINOR_VERSION_1}.${IPE_MINOR_VERSION_2}/ipelets/")
    find_path(IPELET_INSTALL_DIR 
                  NAMES libgoodies.lua goodies.lua
                  PATHS ${INSTALL_PATHS}
                  ENV IPELETPATH
                 )
  else()
    foreach (VER RANGE 28 40)
    string(REPLACE XX ${VER} PATHC "${IPE_LIBRARY_DIR}/ipe/6.0preXX/ipelets/" )
    set(INSTALL_PATHS ${INSTALL_PATHS} ${PATHC})
    endforeach()
    set(INSTALL_PATHS ${INSTALL_PATHS} ${PATHC})
    set(INSTALL_PATHS ${INSTALL_PATHS} /usr/lib64/ipe/6.0/ipelets)
    set(INSTALL_PATHS ${INSTALL_PATHS} /usr/lib/ipe/6.0/ipelets)


    find_library(IPELET_INSTALL_DIR_FILES 
                  NAMES align
                  PATHS ${INSTALL_PATHS}
                        ENV IPELETPATH
                )
    if (IPELET_INSTALL_DIR_FILES)
      get_filename_component(IPELET_INSTALL_DIR ${IPELET_INSTALL_DIR_FILES} PATH)
    endif()                
  endif()
endif()

if(IPE_FOUND)
    message(STATUS "Found Ipe: ${IPE_INCLUDE_DIR} ${IPE_LIBRARY_DIR}")
    if (IPELET_INSTALL_DIR)
      set ( IPELET_INSTALL_DIR ${IPELET_INSTALL_DIR}   CACHE STRING "The folder where ipelets will be installed, relative to CMAKE_INSTALL_PREFIX" )
      message(STATUS "Set Ipelets install dir: ${IPELET_INSTALL_DIR}")
    endif()
endif()
