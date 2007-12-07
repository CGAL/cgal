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

IF(WIN32)
  # In windows, automatic linking is performed, so you do not have to specify the libraries.
  # If you are linking to a dynamic runtime, then you can choose to link to either a static or a
  # dynamic Boost library, the default is to do a static link.  You can alter this for a specific
  # library "whatever" by defining BOOST_WHATEVER_DYN_LINK to force Boost library "whatever" to
  # be linked dynamically.  Alternatively you can force all Boost libraries to dynamic link by
  # defining BOOST_ALL_DYN_LINK.

  # This feature can be disabled for Boost library "whatever" by defining BOOST_WHATEVER_NO_LIB,
  # or for all of Boost by defining BOOST_ALL_NO_LIB.

  # If you want to observe which libraries are being linked against then defining
  # BOOST_LIB_DIAGNOSTIC will cause the auto-linking code to emit a #pragma message each time
  # a library is selected for linking.
  SET(Boost_LIB_DIAGNOSTIC_DEFINITIONS "-DBOOST_LIB_DIAGNOSTIC")
ENDIF()


SET(BOOST_INCLUDE_PATH_DESCRIPTION "directory containing the boost include files. E.g /usr/local/include/boost_1_34_1 or c:\\Program Files\\boost\\boost_1_34_1")

SET(BOOST_DIR_MESSAGE "Set the Boost_INCLUDE_DIR cmake cache entry to the ${BOOST_INCLUDE_PATH_DESCRIPTION}")

SET(BOOST_DIR_SEARCH $ENV{BOOST_ROOT})
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
    "C:/Program Files/boost/boost_1_35_1"
    "C:/Program Files/boost/boost_1_35_0"
    "C:/Program Files/boost/boost_1_34_1"
    "C:/Program Files/boost/boost_1_34_0"
    "C:/Program Files/boost/boost_1_33_1"
    "C:/Program Files/boost/boost_1_33_0"
  )
ENDIF()

#
# Look for an installation.
#
FIND_PATH(Boost_INCLUDE_DIR NAMES boost/config.hpp PATHS ${BOOST_DIR_SEARCH}

  # Help the user find it if we cannot.
  DOC "The ${BOOST_INCLUDE_PATH_DESCRIPTION}"
)

# Assume we didn't find it.
SET(Boost_FOUND 0)

# Now try to get the include and library path.
IF(Boost_INCLUDE_DIR)

  # Compose the boost library path.
  # Note that the user may not have installed any libraries
  # so it is quite possible the Boost_LIBRARY_PATH may not exist.
  SET(Boost_LIBRARY_DIR ${Boost_INCLUDE_DIR})

  # !! ON WHAT PLATFORM CAN THIS WORK !!
  #IF("${Boost_LIBRARY_DIR}" MATCHES "boost-[0-9]+")
  #  GET_FILENAME_COMPONENT(Boost_LIBRARY_DIR ${Boost_LIBRARY_DIR} PATH)
  #ENDIF ("${Boost_LIBRARY_DIR}" MATCHES "boost-[0-9]+")

  # Strip off the trailing "/include" in the path.
  IF("${Boost_LIBRARY_DIR}" MATCHES "/include$")
    GET_FILENAME_COMPONENT(Boost_LIBRARY_DIR ${Boost_LIBRARY_DIR} PATH)
  ENDIF("${Boost_LIBRARY_DIR}" MATCHES "/include$")
  
  IF(EXISTS "${Boost_LIBRARY_DIR}/lib")
    SET (Boost_LIBRARY_DIR ${Boost_LIBRARY_DIR}/lib)
  ELSE()
    IF(EXISTS "${Boost_LIBRARY_DIR}/stage/lib")
      SET(Boost_LIBRARY_DIR ${Boost_LIBRARY_DIR}/stage/lib)
    ELSE()
      SET(Boost_LIBRARY_DIR "")
    ENDIF()
  ENDIF()

  IF(EXISTS "${Boost_INCLUDE_DIR}")
    SET(Boost_INCLUDE_DIRS ${Boost_INCLUDE_DIR})
    # We have found boost. It is possible that the user has not
    # compiled any libraries so we set Boost_FOUND to be true here.
    SET(Boost_FOUND 1)
  ENDIF()

  IF(Boost_LIBRARY_DIR AND EXISTS "${Boost_LIBRARY_DIR}")
    SET(Boost_LIBRARY_DIRS ${Boost_LIBRARY_DIR})
  ENDIF()
ENDIF()

IF(NOT Boost_FOUND)
  IF(NOT Boost2_FIND_QUIETLY)
    MESSAGE(STATUS "Boost was not found. ${BOOST_DIR_MESSAGE}")
  ELSE()
    IF(Boost2_FIND_REQUIRED)
      MESSAGE(FATAL_ERROR "Boost was not found. ${BOOST_DIR_MESSAGE}")
    ENDIF()
  ENDIF()
ENDIF()

