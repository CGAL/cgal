# - Try to find Ipe
# Once done this will define
#
#  IPE_FOUND - system has Ipe
#  IPE_INCLUDE_DIR - the Ipe include directory
#  IPE_LIBRARIES - Link these to use Ipe
#

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

  get_filename_component(IPE_LIBRARY_DIR ${IPE_LIBRARY} PATH)

  
  foreach (VER RANGE 28 40)
  string(REPLACE XX ${VER} PATHC "/usr/lib/ipe/6.0preXX/ipelets/" )
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
               
               

  if(IPE_INCLUDE_DIR AND IPE_LIBRARY)
     set(IPE_FOUND TRUE)
  endif()
endif()

if(IPE_FOUND)
    message(STATUS "Found Ipe: ${IPE_INCLUDE_DIR} ${IPE_LIBRARY_DIR}")
    if (IPELET_INSTALL_DIR)
      set ( IPELET_INSTALL_DIR ${IPELET_INSTALL_DIR}   CACHE STRING "The folder where ipelets will be installed, relative to CMAKE_INSTALL_PREFIX" )
      message(STATUS "Set Ipelets install dir: ${IPELET_INSTALL_DIR}")
    endif()
#~ else()
#~     message(FATAL_ERROR "Could not find Ipe")
endif()
