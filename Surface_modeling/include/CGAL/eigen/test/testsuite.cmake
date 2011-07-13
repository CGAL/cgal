
####################################################################
#
# Usage:
#  - create a new folder, let's call it cdash
#  - in that folder, do:
#    ctest -S path/to/eigen/test/testsuite.cmake[,option1=value1[,option2=value2]]
#
# Options:
#  - EIGEN_CXX: compiler, eg.: g++-4.2
#      default: default c++ compiler
#  - EIGEN_SITE: eg, INRIA-Bdx_pc-gael, or the name of the contributor, etc.
#      default: hostname
#  - EIGEN_BUILD_STRING: a string which identify the system/compiler. It should be formed like that:
#        <OS_name>-<OS_version>-<arch>-<compiler-version>
#      with:
#        <OS_name> = opensuse, debian, osx, windows, cygwin, freebsd, solaris, etc.
#        <OS_version> = 11.1, XP, vista, leopard, etc.
#        <arch> = i386, x86_64, ia64, powerpc, etc.
#        <compiler-version> = gcc-4.3.2, icc-11.0, MSVC-2008, etc.
#  - EIGEN_EXPLICIT_VECTORIZATION: novec, SSE2, Altivec
#       default: SSE2 for x86_64 systems, novec otherwise
#       Its value is automatically appended to EIGEN_BUILD_STRING
#  - EIGEN_CMAKE_DIR: path to cmake executable
#  - EIGEN_MODE: dashboard model, can be Experimental, Nightly, or Continuous
#      default: Nightly
#  - EIGEN_WORK_DIR: directory used to download the source files and make the builds
#      default: folder which contains this script
#  - EIGEN_CMAKE_ARGS: additional arguments passed to cmake
#  - EIGEN_GENERATOR_TYPE: allows to overwrite the generator type
#      default: nmake (windows
#      See http://www.cmake.org/cmake/help/cmake2.6docs.html#section_Generators for a complete
#      list of supported generators.
#  - EIGEN_NO_UPDATE: allows to submit dash boards from local repositories
#      This might be interesting in case you want to submit dashboards
#      including local changes.
#  - CTEST_SOURCE_DIRECTORY: path to eigen's src (use a new and empty folder, not the one you are working on)
#      default: <EIGEN_WORK_DIR>/src
#  - CTEST_BINARY_DIRECTORY: build directory
#      default: <EIGEN_WORK_DIR>/nightly-<EIGEN_CXX>
#
# Here is an example running several compilers on a linux system:
# #!/bin/bash
# ARCH=`uname -m`
# SITE=`hostname`
# VERSION=opensuse-11.1
# WORK_DIR=/home/gael/Coding/eigen/cdash
# # get the last version of the script
# wget http://bitbucket.org/eigen/eigen/raw/tip/test/testsuite.cmake -o $WORK_DIR/testsuite.cmake
# COMMON="ctest -S $WORK_DIR/testsuite.cmake,EIGEN_WORK_DIR=$WORK_DIR,EIGEN_SITE=$SITE,EIGEN_MODE=$1,EIGEN_BUILD_STRING=$OS_VERSION-$ARCH"
# $COMMON-gcc-3.4.6,EIGEN_CXX=g++-3.4
# $COMMON-gcc-4.0.1,EIGEN_CXX=g++-4.0.1
# $COMMON-gcc-4.3.2,EIGEN_CXX=g++-4.3,EIGEN_EXPLICIT_VECTORIZATION=novec
# $COMMON-gcc-4.3.2,EIGEN_CXX=g++-4.3,EIGEN_EXPLICIT_VECTORIZATION=SSE2
# $COMMON-icc-11.0,EIGEN_CXX=icpc
#
####################################################################

# process the arguments

set(ARGLIST ${CTEST_SCRIPT_ARG})
while(${ARGLIST} MATCHES  ".+.*")

  # pick first
  string(REGEX MATCH "([^,]*)(,.*)?" DUMMY ${ARGLIST})
  SET(TOP ${CMAKE_MATCH_1})

  # remove first
  string(REGEX MATCHALL "[^,]*,(.*)" DUMMY ${ARGLIST})
  SET(ARGLIST ${CMAKE_MATCH_1})

  # decompose as a pair key=value
  string(REGEX MATCH "([^=]*)(=.*)?" DUMMY ${TOP})
  SET(KEY ${CMAKE_MATCH_1})

  string(REGEX MATCH "[^=]*=(.*)" DUMMY ${TOP})
  SET(VALUE ${CMAKE_MATCH_1})

  # set the variable to the specified value
  if(VALUE)
    SET(${KEY} ${VALUE})
  else(VALUE)
    SET(${KEY} ON)
  endif(VALUE)

endwhile(${ARGLIST} MATCHES ".+.*")

####################################################################
# Automatically set some user variables if they have not been defined manually
####################################################################
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)

if(NOT EIGEN_SITE)
  site_name(EIGEN_SITE)
endif(NOT EIGEN_SITE)

if(NOT EIGEN_CMAKE_DIR)
  SET(EIGEN_CMAKE_DIR "")
endif(NOT EIGEN_CMAKE_DIR)

if(NOT EIGEN_BUILD_STRING)

  # let's try to find all information we need to make the build string ourself

  # OS
  build_name(EIGEN_OS_VERSION)

  # arch
  set(EIGEN_ARCH ${CMAKE_SYSTEM_PROCESSOR})
  if(WIN32)
    set(EIGEN_ARCH $ENV{PROCESSOR_ARCHITECTURE})
  else(WIN32)
    execute_process(COMMAND uname -m OUTPUT_VARIABLE EIGEN_ARCH OUTPUT_STRIP_TRAILING_WHITESPACE)
  endif(WIN32)

  set(EIGEN_BUILD_STRING ${EIGEN_OS_VERSION}${EIGEN_ARCH}-${EIGEN_CXX})

endif(NOT EIGEN_BUILD_STRING)

if(DEFINED EIGEN_EXPLICIT_VECTORIZATION)
  set(EIGEN_BUILD_STRING ${EIGEN_BUILD_STRING}-${EIGEN_EXPLICIT_VECTORIZATION})
endif(DEFINED EIGEN_EXPLICIT_VECTORIZATION)

if(NOT EIGEN_WORK_DIR)
  set(EIGEN_WORK_DIR ${CTEST_SCRIPT_DIRECTORY})
endif(NOT EIGEN_WORK_DIR)

if(NOT CTEST_SOURCE_DIRECTORY)
  SET (CTEST_SOURCE_DIRECTORY "${EIGEN_WORK_DIR}/src")
endif(NOT CTEST_SOURCE_DIRECTORY)

if(NOT CTEST_BINARY_DIRECTORY)
  SET (CTEST_BINARY_DIRECTORY "${EIGEN_WORK_DIR}/nightly_${EIGEN_CXX}")
endif(NOT CTEST_BINARY_DIRECTORY)

if(NOT EIGEN_MODE)
  set(EIGEN_MODE Nightly)
endif(NOT EIGEN_MODE)

## mandatory variables (the default should be ok in most cases):

if(NOT EIGEN_NO_UPDATE)
  SET (CTEST_CVS_COMMAND "hg")
  SET (CTEST_CVS_CHECKOUT "${CTEST_CVS_COMMAND} clone http://bitbucket.org/eigen/eigen \"${CTEST_SOURCE_DIRECTORY}\"")
  SET(CTEST_BACKUP_AND_RESTORE TRUE) # the backup is CVS related ...
endif(NOT EIGEN_NO_UPDATE)

# which ctest command to use for running the dashboard
SET (CTEST_COMMAND "${EIGEN_CMAKE_DIR}ctest -D ${EIGEN_MODE} --no-compress-output")
if($ENV{EIGEN_CTEST_ARGS})
SET (CTEST_COMMAND "${CTEST_COMMAND} $ENV{EIGEN_CTEST_ARGS}")
endif($ENV{EIGEN_CTEST_ARGS})
# what cmake command to use for configuring this dashboard
SET (CTEST_CMAKE_COMMAND "${EIGEN_CMAKE_DIR}cmake -DEIGEN_LEAVE_TEST_IN_ALL_TARGET=ON")

####################################################################
# The values in this section are optional you can either
# have them or leave them commented out
####################################################################

# this make sure we get consistent outputs
SET($ENV{LC_MESSAGES} "en_EN")

# should ctest wipe the binary tree before running
SET(CTEST_START_WITH_EMPTY_BINARY_DIRECTORY TRUE)

# raise the warning/error limit
set(CTEST_CUSTOM_MAXIMUM_NUMBER_OF_WARNINGS "33331")
set(CTEST_CUSTOM_MAXIMUM_NUMBER_OF_ERRORS "33331")

# this is the initial cache to use for the binary tree, be careful to escape
# any quotes inside of this string if you use it
if(WIN32 AND NOT UNIX)
  #message(SEND_ERROR "win32")
  if(EIGEN_GENERATOR_TYPE)
    set(CTEST_CMAKE_COMMAND "${CTEST_CMAKE_COMMAND} -G \"${EIGEN_GENERATOR_TYPE}\"")
    SET (CTEST_INITIAL_CACHE "
      CMAKE_BUILD_TYPE:STRING=Release
      BUILDNAME:STRING=${EIGEN_BUILD_STRING}
      SITE:STRING=${EIGEN_SITE}
    ")
  else(EIGEN_GENERATOR_TYPE)
    set(CTEST_CMAKE_COMMAND "${CTEST_CMAKE_COMMAND} -G \"NMake Makefiles\" -DCMAKE_MAKE_PROGRAM=nmake")
    SET (CTEST_INITIAL_CACHE "
      MAKECOMMAND:STRING=nmake /i
      CMAKE_MAKE_PROGRAM:FILEPATH=nmake
      CMAKE_GENERATOR:INTERNAL=NMake Makefiles
      CMAKE_BUILD_TYPE:STRING=Release
      BUILDNAME:STRING=${EIGEN_BUILD_STRING}
      SITE:STRING=${EIGEN_SITE}
    ")
  endif(EIGEN_GENERATOR_TYPE)
else(WIN32 AND NOT UNIX)
  SET (CTEST_INITIAL_CACHE "
    BUILDNAME:STRING=${EIGEN_BUILD_STRING}
    SITE:STRING=${EIGEN_SITE}
  ")
endif(WIN32 AND NOT UNIX)

# set any extra environment variables to use during the execution of the script here:
# setting this variable on windows machines causes trouble ...

if(EIGEN_CXX AND NOT WIN32)
  set(CTEST_ENVIRONMENT "CXX=${EIGEN_CXX}")
endif(EIGEN_CXX AND NOT WIN32)

if(DEFINED EIGEN_EXPLICIT_VECTORIZATION)
  if(EIGEN_EXPLICIT_VECTORIZATION MATCHES SSE2)
    set(CTEST_CMAKE_COMMAND "${CTEST_CMAKE_COMMAND} -DEIGEN_TEST_SSE2=ON")
  elseif(EIGEN_EXPLICIT_VECTORIZATION MATCHES SSE3)
    set(CTEST_CMAKE_COMMAND "${CTEST_CMAKE_COMMAND} -DEIGEN_TEST_SSE2=ON -DEIGEN_TEST_SSE3=ON")
  elseif(EIGEN_EXPLICIT_VECTORIZATION MATCHES SSSE3)
    set(CTEST_CMAKE_COMMAND "${CTEST_CMAKE_COMMAND} -DEIGEN_TEST_SSE2=ON -DEIGEN_TEST_SSE3=ON -DEIGEN_TEST_SSSE3=ON")  
  elseif(EIGEN_EXPLICIT_VECTORIZATION MATCHES SSE4_1)
    set(CTEST_CMAKE_COMMAND "${CTEST_CMAKE_COMMAND} -DEIGEN_TEST_SSE2=ON -DEIGEN_TEST_SSE3=ON -DEIGEN_TEST_SSSE3=ON -DEIGEN_TEST_SSE4_1=ON")
  elseif(EIGEN_EXPLICIT_VECTORIZATION MATCHES SSE4_2)
    set(CTEST_CMAKE_COMMAND "${CTEST_CMAKE_COMMAND} -DEIGEN_TEST_SSE2=ON -DEIGEN_TEST_SSE3=ON -DEIGEN_TEST_SSSE3=ON -DEIGEN_TEST_SSE4_1=ON -DEIGEN_TEST_SSE4_2=ON")
  elseif(EIGEN_EXPLICIT_VECTORIZATION MATCHES Altivec)
    set(CTEST_CMAKE_COMMAND "${CTEST_CMAKE_COMMAND} -DEIGEN_TEST_ALTIVEC=ON")
  elseif(EIGEN_EXPLICIT_VECTORIZATION MATCHES novec)
    set(CTEST_CMAKE_COMMAND "${CTEST_CMAKE_COMMAND} -DEIGEN_TEST_NO_EXPLICIT_VECTORIZATION=ON")
  else(EIGEN_EXPLICIT_VECTORIZATION MATCHES SSE2)
    message(FATAL_ERROR "Invalid value for EIGEN_EXPLICIT_VECTORIZATION (${EIGEN_EXPLICIT_VECTORIZATION}), must be: novec, SSE2, SSE3, Altivec")
  endif(EIGEN_EXPLICIT_VECTORIZATION MATCHES SSE2)
endif(DEFINED EIGEN_EXPLICIT_VECTORIZATION)

if(DEFINED EIGEN_CMAKE_ARGS)
  set(CTEST_CMAKE_COMMAND "${CTEST_CMAKE_COMMAND} ${EIGEN_CMAKE_ARGS}")
endif(DEFINED EIGEN_CMAKE_ARGS)
