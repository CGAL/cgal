# Copyright (c) 2012
# Utrecht University (The Netherlands),
# ETH Zurich (Switzerland),
# INRIA Sophia-Antipolis (France),
# Max-Planck-Institute Saarbruecken (Germany),
# and Tel-Aviv University (Israel).  All rights reserved.
#
# This file is part of CGAL (www.cgal.org); you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public License as
# published by the Free Software Foundation; version 3 of the License,
# or (at your option) any later version.
#
# Licensees holding a valid commercial license may use this file in
# accordance with the commercial license agreement provided with the software.
#
# This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
# WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
#
# $URL: svn+ssh://fcacciola@scm.gforge.inria.fr/svn/cgal/trunk/Scripts/scripts/cgal_create_makefile $
# $Id: cgal_create_makefile 36976 2007-03-09 22:53:24Z reichel $
#
# Author(s)     : Eric Berberich <eric@mpi-inf.mpg.de>

cmake_minimum_required(VERSION 2.6.2)

message(STATUS "Create CMakeLists.txt")

# message(STATUS "Repeat command line options: ${OPTIONS}")

set(PROJECT CGAL) #`basename $PWD`
set(SINGLE_SOURCE "")
set(CGAL_COMPONENTS "QT4 GMP MPFR RS3")
set(WITH_QT3 FALSE)
set(WITH_QT4 FALSE)
set(WITH_ALL_PRECONFIGURED_LIBS TRUE)
set(BOOST_COMPONENTS "thread")

# TODO enable_testing()?


### Delete file if it exists

if (EXISTS CMakeLists.txt)
  file(RENAME CMakeLists.txt CMakeLists.bak)
endif()

### Parse options

# TODO

#-s
#-c/-p
#-b
#-t? # for testing

### Start to write file

file(APPEND CMakeLists.txt "# Created by the script cgal_generate_cmake_script\n" )
file(APPEND CMakeLists.txt "# This is the CMake script for compiling a set of CGAL applications.\n\n" )

if ( "${SINGLE_SOURCE}" STREQUAL "" ) 
  file(APPEND CMakeLists.txt "project(${PROJECT})\n\n")
else()
  file(APPEND CMakeLists.txt "project(${SINGLE_SOURCE})\n\n")
endif()


file(APPEND CMakeLists.txt
"cmake_minimum_required(VERSION 2.6.2)
if(\"\${CMAKE_MAJOR_VERSION}.\${CMAKE_MINOR_VERSION}\" VERSION_GREATER 2.6)
  if(\"\${CMAKE_MAJOR_VERSION}.\${CMAKE_MINOR_VERSION}.\${CMAKE_PATCH_VERSION}\" VERSION_GREATER 2.8.3)
    cmake_policy(VERSION 2.8.4)
  else()
    cmake_policy(VERSION 2.6)
  endif()
endif()

set( CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true )
 
if ( COMMAND cmake_policy )

  cmake_policy( SET CMP0003 NEW )  

endif()

# CGAL and its components
")

# TODO testing?
#if [ -n "$ENABLE_CTEST" ]; then
#  echo "enable_testing()"
#fi

# TODO case of components, 

foreach( component ${CGAL_COMPONENTS})
  # detect qt3, qt4, and ALL_PRECONFIGURED_LIBS
  if ( ${component} STREQUAL "Qt3" )
     set(WITH_QT3 TRUE)
  endif()
  if ( ${component} STREQUAL "Qt4" )
     set(WITH_QT4 TRUE)
  endif()
endforeach()

# TODO replace " "  with ":"

if ( WITH_ALL_PRECONFIGURED_LIBS )
  set(CGAL_COMPONENTS "${CGAL_COMPONENTS}:ALL_PRECONFIGURED_LIBS")
endif()

file(APPEND CMakeLists.txt "find_package( CGAL QUIET COMPONENTS ${CGAL_COMPONENTS} )\n\n")

file(APPEND CMakeLists.txt 
"if ( NOT CGAL_FOUND )

  message(STATUS \"This project requires the CGAL library, and will not be compiled.\")
  return()  

endif()

# include helper file
include( \${CGAL_USE_FILE} )

# Boost and its components
")


### Boost and its components

# TODO replace " "  with ":"

if ( ${BOOST_COMPONENTS} STREQUAL "")
  file(APPEND CMakeLists.txt "find_package( Boost REQUIRED )\n\n")
else()
  file(APPEND CMakeLists.txt "find_package( Boost REQUIRED COMPONENTS ${BOOST_COMPONENTS} )\n\n")
endif()

file(APPEND CMakeLists.txt 
"if ( NOT Boost_FOUND )

  message(STATUS \"This project requires the Boost library, and will not be compiled.\")

  return()  

endif()

")

if ( EXISTS include ) 
  file(APPEND CMakeLists.txt
"# include for local directory
include_directories( BEFORE include )\n\n")
endif()

if ( EXISTS ../../include ) 
  file(APPEND CMakeLists.txt
"# includes for local package
include_directories( BEFORE ../../include )\n\n")
endif()

if ( EXISTS ../include ) 
  file(APPEND CMakeLists.txt
"# includes for local package
include_directories( BEFORE ../include )\n\n")
endif()
