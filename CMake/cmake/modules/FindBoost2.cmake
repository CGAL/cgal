#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# This is FindBoost2, a fork of FindBoost which works also with the default directory structure on windows.
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# - Find the Boost includes and libraries.
# The following variables are set if Boost is found.  If Boost is not
# found, Boost_FOUND is set to false.
#  Boost_FOUND        - True when the Boost include directory is found.
#  Boost_INCLUDE_DIRS - the path to where the boost include files are.
#  Boost_LIBRARY_DIRS - The path to where the boost library files are.
#  Boost_LIB_DIAGNOSTIC_DEFINITIONS - Only set if using Windows.

# ----------------------------------------------------------------------------
# If you have installed Boost in a non-standard location or you have
# just staged the boost files using bjam then you have three
# options. In the following comments, it is assumed that <Your Path>
# points to the root directory of the include directory of Boost. e.g
# If you have put boost in C:\development\Boost then <Your Path> is
# "C:/development/Boost" and in this directory there will be two
# directories called "include" and "lib".
# 1) After CMake runs, set Boost_INCLUDE_DIR to <Your Path>/include/boost<-version>
# 2) Use CMAKE_INCLUDE_PATH to set a path to <Your Path>/include. This will allow FIND_PATH()
#    to locate Boost_INCLUDE_DIR by utilizing the PATH_SUFFIXES option. e.g.
#    SET(CMAKE_INCLUDE_PATH ${CMAKE_INCLUDE_PATH} "<Your Path>/include")
# 3) Set an environment variable called ${BOOSTROOT} that points to the root of where you have
#    installed Boost, e.g. <Your Path>. It is assumed that there is at least a subdirectory called
#    include in this path.
#
# Note:
#  1) If you are just using the boost headers, then you do not need to use
#     Boost_LIBRARY_DIRS in your CMakeLists.txt file.
#  2) If Boost has not been installed, then when setting Boost_LIBRARY_DIRS
#     the script will look for /lib first and, if this fails, then for /stage/lib.
#
# Usage:
# In your CMakeLists.txt file do something like this:
# ...
# # Boost
# FIND_PACKAGE(Boost)
# ...
# INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})
# LINK_DIRECTORIES(${Boost_LIBRARY_DIRS})
#
# In Windows, we make the assumption that, if the Boost files are installed, the default directory
# will be C:\boost.

#
# TODO:
#
# 1) Automatically find the Boost library files and eliminate the need
#    to use Link Directories.
#

SET(BOOST_INCLUDE_PATH_DESCRIPTION "directory containing the boost include files. E.g /usr/local/include/boost_1_34_1 or c:\\Program Files\\boost\\boost_1_34_1")

SET(BOOST_LIB_PATH_DESCRIPTION "directory containing the boost library files. E.g /usr/local/lib/boost_1_34_1 or c:\\Program Files\\boost\\boost_1_34_1\\lib")

SET(BOOST_DIR_MESSAGE "Set the Boost_INCLUDE_DIR cmake cache entry to the directory containing the boost include files.")

SET(BOOST_DIR_SEARCH $ENV{BOOST_ROOT})
if ( NOT BOOST_DIR_SEARCH ) 
  SET(BOOST_DIR_SEARCH $ENV{BOOSTROOT} ) 
endif()

IF(BOOST_DIR_SEARCH)
  FILE(TO_CMAKE_PATH ${BOOST_DIR_SEARCH} BOOST_DIR_SEARCH)
  
  # In Windows, BOOST_ROOT refers directory to the folder having boost/<headers>.hpp
  # but on Linux, BOOST_ROOT typically refers to a folder having incude/boost/<headers>.hpp
  SET(BOOST_DIR_SEARCH ${BOOST_DIR_SEARCH} {BOOST_DIR_SEARCH}/include)
ENDIF()

# The windows installer uses C:\\Program Files\\boost\\boost_XX_YY_ZZ by default
IF(WIN32)
  SET(BOOST_DIR_SEARCH
    ${BOOST_DIR_SEARCH}
    
    "$ENV{ProgramFiles}/boost/boost_1_35_1"
    "$ENV{ProgramFiles}/boost/boost_1_35_0"
    "$ENV{ProgramFiles}/boost/boost_1_34_1"
    "$ENV{ProgramFiles}/boost/boost_1_34_0"
    "$ENV{ProgramFiles}/boost/boost_1_33_1"
    "$ENV{ProgramFiles}/boost/boost_1_33_0"
    
    "$ENV{ProgramFiles}/boost_1_35_1"
    "$ENV{ProgramFiles}/boost_1_35_0"
    "$ENV{ProgramFiles}/boost_1_34_1"
    "$ENV{ProgramFiles}/boost_1_34_0"
    "$ENV{ProgramFiles}/boost_1_33_1"
    "$ENV{ProgramFiles}/boost_1_33_0"
  )
ENDIF()

#
# Look for an installation.
#
FIND_PATH(Boost_INCLUDE_DIRS NAMES boost/config.hpp PATHS ${BOOST_DIR_SEARCH} DOC "The ${BOOST_INCLUDE_PATH_DESCRIPTION}" )

# Now try to get the include and library path.
IF(Boost_INCLUDE_DIRS)
  
  # Compose the boost library path.
  # Note that the user may not have installed any libraries
  # so it is quite possible the Boost_LIBRARY_PATH may not exist.
  SET(Boost_LIBRARY_DIR ${Boost_INCLUDE_DIRS} )

  # !! ON WHAT PLATFORM CAN THIS WORK !!
  #IF("${Boost_LIBRARY_DIR}" MATCHES "boost-[0-9]+")
  #  GET_FILENAME_COMPONENT(Boost_LIBRARY_DIR ${Boost_LIBRARY_DIR} PATH)
  #ENDIF ("${Boost_LIBRARY_DIR}" MATCHES "boost-[0-9]+")

  # Strip off the trailing "/include" in the path.
  IF("${Boost_LIBRARY_DIR}" MATCHES "/include$")
    GET_FILENAME_COMPONENT(Boost_LIBRARY_DIR ${Boost_LIBRARY_DIR} PATH)
  ENDIF()
  
  IF(EXISTS "${Boost_LIBRARY_DIR}/lib")
    SET (Boost_LIBRARY_DIR ${Boost_LIBRARY_DIR}/lib)
  ELSE()
    IF(EXISTS "${Boost_LIBRARY_DIR}/stage/lib")
      SET(Boost_LIBRARY_DIR ${Boost_LIBRARY_DIR}/stage/lib)
    ELSE()
      SET(Boost_LIBRARY_DIR "")
    ENDIF()
  ENDIF()

  IF(Boost_LIBRARY_DIR AND EXISTS "${Boost_LIBRARY_DIR}")
    SET(Boost_LIBRARY_DIRS ${Boost_LIBRARY_DIR} CACHE PATH "The ${BOOST_LIB_PATH_DESCRIPTION}" )
  ENDIF()
ENDIF()

INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(Boost "Boost was not found. ${BOOST_DIR_MESSAGE}" Boost_INCLUDE_DIRS )
SET(Boost_FOUND ${BOOST_FOUND})

MARK_AS_ADVANCED(Boost_INCLUDE_DIRS)
MARK_AS_ADVANCED(Boost_LIBRARY_DIRS)
