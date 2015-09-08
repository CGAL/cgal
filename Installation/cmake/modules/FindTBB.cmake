# Locate Intel Threading Building Blocks include paths and libraries
# FindTBB.cmake can be found at https://code.google.com/p/findtbb/
# Written by Hannes Hofmann <hannes.hofmann _at_ informatik.uni-erlangen.de>
# Improvements by Gino van den Bergen <gino _at_ dtecta.com>,
#   Florian Uhlig <F.Uhlig _at_ gsi.de>,
#   Jiri Marsik <jiri.marsik89 _at_ gmail.com>

# The MIT License
#
# Copyright (c) 2011 Hannes Hofmann
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.

# GvdB: This module uses the environment variable TBB_ARCH_PLATFORM which defines architecture and compiler.
#   e.g. "ia32/vc8" or "em64t/cc4.1.0_libc2.4_kernel2.6.16.21"
#   TBB_ARCH_PLATFORM is set by the build script tbbvars[.bat|.sh|.csh], which can be found
#   in the TBB installation directory (TBB_INSTALL_DIR).
#
# GvdB: Mac OS X distribution places libraries directly in lib directory.
#
# For backwards compatibility, you may explicitely set the CMake variables TBB_ARCHITECTURE and TBB_COMPILER.
# TBB_ARCHITECTURE [ ia32 | em64t | itanium ]
#   which architecture to use
# TBB_COMPILER e.g. vc9 or cc3.2.3_libc2.3.2_kernel2.4.21 or cc4.0.1_os10.4.9
#   which compiler to use (detected automatically on Windows)

# This module respects
# TBB_INSTALL_DIR or $ENV{TBB21_INSTALL_DIR} or $ENV{TBB_INSTALL_DIR}

# This module defines
# TBB_INCLUDE_DIRS, where to find task_scheduler_init.h, etc.
# TBB_LIBRARY_DIRS, where to find TBB libraries
# TBB_INSTALL_DIR, the base TBB install directory.
# TBB_LIBRARIES, all the following TBB libraries (both release and debug versions, using "optimized" and "debug" CMake keywords). Note that if the debug versions are not found, the release versions will be used instead for the debug mode.
#   TBB_RELEASE_LIBRARY, the TBB release library
#   TBB_MALLOC_RELEASE_LIBRARY, the TBB release malloc library
#   TBB_DEBUG_LIBRARY, the TBB debug library
#   TBB_MALLOC_DEBUG_LIBRARY, the TBB debug malloc library
# TBB_FOUND, If false, don't try to use TBB.
# TBB_INTERFACE_VERSION, as defined in tbb/tbb_stddef.h
# TBB_MALLOCPROXY_DEBUG_LIBRARY, the TBB debug malloc_proxy library (not included in TBB_LIBRARIES since it's optionnal)
# TBB_MALLOCPROXY_RELEASE_LIBRARY, the TBB release malloc_proxy library (not included in TBB_LIBRARIES since it's optionnal)

include(CheckCXXSourceCompiles)

# Usage:
#   try_TBB_with_pthread(<result_var_name> [additional linker args...])
function(try_TBB_with_pthread result_var)
    set(TBB_try_ts_source "
          #include <tbb/enumerable_thread_specific.h>
          int main() {
            tbb::enumerable_thread_specific<
              bool*,
              tbb::cache_aligned_allocator<bool*>,
              tbb::ets_key_per_instance> grid;
          }
          ")
    set(CMAKE_REQUIRED_LIBRARIES ${ALL_TBB_LIBRARIES} ${ARGN})
    set(CMAKE_REQUIRED_INCLUDES ${TBB_INCLUDE_DIR})
    check_cxx_source_compiles("${TBB_try_ts_source}" ${result_var})
    set(${result_var} ${${result_var}} PARENT_SCOPE)
endfunction(try_TBB_with_pthread)

if (WIN32)
    # has em64t/vc8 em64t/vc9
    # has ia32/vc7.1 ia32/vc8 ia32/vc9
    set(_TBB_DEFAULT_INSTALL_DIR "C:/Program Files/Intel/TBB" "C:/Program Files (x86)/Intel/TBB")
    set(_TBB_LIB_RELEASE_NAME "tbb")
    set(_TBB_LIB_MALLOC_RELEASE_NAME "${_TBB_LIB_RELEASE_NAME}malloc")
    set(_TBB_LIB_MALLOCPROXY_RELEASE_NAME "${_TBB_LIB_RELEASE_NAME}malloc_proxy")
    set(_TBB_LIB_DEBUG_NAME "${_TBB_LIB_RELEASE_NAME}_debug")
    set(_TBB_LIB_MALLOC_DEBUG_NAME "${_TBB_LIB_MALLOC_RELEASE_NAME}_debug")
    set(_TBB_LIB_MALLOCPROXY_DEBUG_NAME "${_TBB_LIB_MALLOCPROXY_RELEASE_NAME}_debug")
    if (MSVC71)
        set (_TBB_COMPILER "vc7.1")
    endif(MSVC71)
    if (MSVC80)
        set(_TBB_COMPILER "vc8")
    endif(MSVC80)
    if (MSVC90)
        set(_TBB_COMPILER "vc9")
    endif(MSVC90)
    if(MSVC10)
        set(_TBB_COMPILER "vc10")
    endif(MSVC10)
    if(MSVC11)
        set(_TBB_COMPILER "vc11")
    endif(MSVC11)
    if(MSVC12)
	set(_TBB_COMPILER "vc12")
    endif(MSVC12)
    #note there was no MSVC13
    if(MSVC14)
	if(RUNNING_CGAL_AUTO_TEST)
	    set (TBB_FOUND "NO")
	    return()#binaries for TBB not publicly available when CGAL-4.7 is published
	endif(RUNNING_CGAL_AUTO_TEST)
	message(STATUS "[Warning] FindTBB.cmake: TBB 4.4 (latest available when CGAL-4.7 is published) does not provide support for MSVC 2015.")
    endif(MSVC14)
    # Todo: add other Windows compilers such as ICL.
    set(_TBB_ARCHITECTURE ${TBB_ARCHITECTURE})
endif (WIN32)

if (UNIX)
    if (APPLE)
        # MAC
        set(_TBB_DEFAULT_INSTALL_DIR "/Library/Frameworks/Intel_TBB.framework/Versions")
        # libs: libtbb.dylib, libtbbmalloc.dylib, *_debug
        set(_TBB_LIB_RELEASE_NAME "tbb")
        set(_TBB_LIB_MALLOC_RELEASE_NAME "${_TBB_LIB_RELEASE_NAME}malloc")
        #set(_TBB_LIB_MALLOCPROXY_RELEASE_NAME "${_TBB_LIB_RELEASE_NAME}malloc_proxy")
        set(_TBB_LIB_DEBUG_NAME "${_TBB_LIB_RELEASE_NAME}_debug")
        set(_TBB_LIB_MALLOC_DEBUG_NAME "${_TBB_LIB_MALLOC_RELEASE_NAME}_debug")
        #set(_TBB_LIB_MALLOCPROXY_DEBUG_NAME "${_TBB_LIB_MALLOCPROXY_RELEASE_NAME}_debug")
        # default flavor on apple: ia32/cc4.0.1_os10.4.9
        # Jiri: There is no reason to presume there is only one flavor and
        #       that user's setting of variables should be ignored.
        if(NOT TBB_COMPILER)
            set(_TBB_COMPILER "cc4.0.1_os10.4.9")
        elseif (NOT TBB_COMPILER)
            set(_TBB_COMPILER ${TBB_COMPILER})
        endif(NOT TBB_COMPILER)
        if(NOT TBB_ARCHITECTURE)
            set(_TBB_ARCHITECTURE "ia32")
        elseif(NOT TBB_ARCHITECTURE)
            set(_TBB_ARCHITECTURE ${TBB_ARCHITECTURE})
        endif(NOT TBB_ARCHITECTURE)
    else (APPLE)
        # LINUX
        set(_TBB_DEFAULT_INSTALL_DIR "/opt/intel/tbb" "/usr/local/include" "/usr/include")
        set(_TBB_LIB_RELEASE_NAME "tbb")
        set(_TBB_LIB_MALLOC_RELEASE_NAME "${_TBB_LIB_RELEASE_NAME}malloc")
        set(_TBB_LIB_MALLOCPROXY_RELEASE_NAME "${_TBB_LIB_RELEASE_NAME}malloc_proxy")
        set(_TBB_LIB_DEBUG_NAME "${_TBB_LIB_RELEASE_NAME}_debug")
        set(_TBB_LIB_MALLOC_DEBUG_NAME "${_TBB_LIB_MALLOC_RELEASE_NAME}_debug")
        set(_TBB_LIB_MALLOCPROXY_DEBUG_NAME "${_TBB_LIB_MALLOCPROXY_RELEASE_NAME}_debug")
        # has em64t/cc3.2.3_libc2.3.2_kernel2.4.21 em64t/cc3.3.3_libc2.3.3_kernel2.6.5 em64t/cc3.4.3_libc2.3.4_kernel2.6.9 em64t/cc4.1.0_libc2.4_kernel2.6.16.21
        # has ia32/*
        # has itanium/*
        set(_TBB_COMPILER ${TBB_COMPILER})
        set(_TBB_ARCHITECTURE ${TBB_ARCHITECTURE})
    endif (APPLE)
endif (UNIX)

if (CMAKE_SYSTEM MATCHES "SunOS.*")
# SUN
# not yet supported
# has em64t/cc3.4.3_kernel5.10
# has ia32/*
endif (CMAKE_SYSTEM MATCHES "SunOS.*")


#-- Clear the public variables
set (TBB_FOUND "NO")


#-- Find TBB install dir and set ${_TBB_INSTALL_DIR} and cached ${TBB_INSTALL_DIR}
# first: use CMake variable TBB_INSTALL_DIR
if (TBB_INSTALL_DIR)
    set (_TBB_INSTALL_DIR ${TBB_INSTALL_DIR})
endif (TBB_INSTALL_DIR)
# second: use environment variable
if (NOT _TBB_INSTALL_DIR)
    if (NOT "$ENV{TBBROOT}" STREQUAL "")
        set (_TBB_INSTALL_DIR $ENV{TBBROOT})
    endif (NOT "$ENV{TBBROOT}" STREQUAL "")
    if (NOT "$ENV{TBB_INSTALL_DIR}" STREQUAL "")
        set (_TBB_INSTALL_DIR $ENV{TBB_INSTALL_DIR})
    endif (NOT "$ENV{TBB_INSTALL_DIR}" STREQUAL "")
    # Intel recommends setting TBB21_INSTALL_DIR
    if (NOT "$ENV{TBB21_INSTALL_DIR}" STREQUAL "")
        set (_TBB_INSTALL_DIR $ENV{TBB21_INSTALL_DIR})
    endif (NOT "$ENV{TBB21_INSTALL_DIR}" STREQUAL "")
    if (NOT "$ENV{TBB22_INSTALL_DIR}" STREQUAL "")
        set (_TBB_INSTALL_DIR $ENV{TBB22_INSTALL_DIR})
    endif (NOT "$ENV{TBB22_INSTALL_DIR}" STREQUAL "")
    if (NOT "$ENV{TBB30_INSTALL_DIR}" STREQUAL "")
        set (_TBB_INSTALL_DIR $ENV{TBB30_INSTALL_DIR})
    endif (NOT "$ENV{TBB30_INSTALL_DIR}" STREQUAL "")
endif (NOT _TBB_INSTALL_DIR)
# third: try to find path automatically
if (NOT _TBB_INSTALL_DIR)
    if (_TBB_DEFAULT_INSTALL_DIR)
        set (_TBB_INSTALL_DIR ${_TBB_DEFAULT_INSTALL_DIR})
    endif (_TBB_DEFAULT_INSTALL_DIR)
endif (NOT _TBB_INSTALL_DIR)
# sanity check
if (NOT _TBB_INSTALL_DIR)
    message ("ERROR: Unable to find Intel TBB install directory. ${_TBB_INSTALL_DIR}")
else (NOT _TBB_INSTALL_DIR)
# finally: set the cached CMake variable TBB_INSTALL_DIR
if (NOT TBB_INSTALL_DIR)
    set (TBB_INSTALL_DIR ${_TBB_INSTALL_DIR} CACHE PATH "Intel TBB install directory")
    mark_as_advanced(TBB_INSTALL_DIR)
endif (NOT TBB_INSTALL_DIR)


#-- A macro to rewrite the paths of the library. This is necessary, because
#   find_library() always found the em64t/vc9 version of the TBB libs
macro(TBB_CORRECT_LIB_DIR var_name)
#    if (NOT "${_TBB_ARCHITECTURE}" STREQUAL "em64t")
        string(REPLACE em64t "${_TBB_ARCHITECTURE}" ${var_name} ${${var_name}})
#    endif (NOT "${_TBB_ARCHITECTURE}" STREQUAL "em64t")
    string(REPLACE ia32 "${_TBB_ARCHITECTURE}" ${var_name} ${${var_name}})
    string(REPLACE vc7.1 "${_TBB_COMPILER}" ${var_name} ${${var_name}})
    string(REPLACE vc8 "${_TBB_COMPILER}" ${var_name} ${${var_name}})
    string(REPLACE vc9 "${_TBB_COMPILER}" ${var_name} ${${var_name}})
    string(REPLACE vc10 "${_TBB_COMPILER}" ${var_name} ${${var_name}})
    string(REPLACE vc11 "${_TBB_COMPILER}" ${var_name} ${${var_name}})
endmacro(TBB_CORRECT_LIB_DIR var_content)


#-- Look for include directory and set ${TBB_INCLUDE_DIR}
set (TBB_INC_SEARCH_DIR ${_TBB_INSTALL_DIR}/include)
# Jiri: tbbvars now sets the CPATH environment variable to the directory
#       containing the headers.
#  LR: search first with NO_DEFAULT_PATH...
find_path(TBB_INCLUDE_DIR
    tbb/task_scheduler_init.h
    PATHS ${TBB_INC_SEARCH_DIR} ENV CPATH
    NO_DEFAULT_PATH
)
if(NOT TBB_INCLUDE_DIR)
#  LR: ... and then search again with NO_DEFAULT_PATH if nothing was found in
#  hinted paths
  find_path(TBB_INCLUDE_DIR
      tbb/task_scheduler_init.h
      PATHS ${TBB_INC_SEARCH_DIR} ENV CPATH
  )
endif()
mark_as_advanced(TBB_INCLUDE_DIR)


#-- Look for libraries
# GvdB: $ENV{TBB_ARCH_PLATFORM} is set by the build script tbbvars[.bat|.sh|.csh]
if (NOT $ENV{TBB_ARCH_PLATFORM} STREQUAL "")
    set (_TBB_LIBRARY_DIR
         ${_TBB_INSTALL_DIR}/lib/$ENV{TBB_ARCH_PLATFORM}
         ${_TBB_INSTALL_DIR}/$ENV{TBB_ARCH_PLATFORM}/lib
        )
endif (NOT $ENV{TBB_ARCH_PLATFORM} STREQUAL "")
# Jiri: This block isn't mutually exclusive with the previous one
#       (hence no else), instead I test if the user really specified
#       the variables in question.
if ((NOT ${TBB_ARCHITECTURE} STREQUAL "") AND (NOT ${TBB_COMPILER} STREQUAL ""))
    # HH: deprecated
    message(STATUS "[Warning] FindTBB.cmake: The use of TBB_ARCHITECTURE and TBB_COMPILER is deprecated and may not be supported in future versions. Please set \$ENV{TBB_ARCH_PLATFORM} (using tbbvars.[bat|csh|sh]).")
    # Jiri: It doesn't hurt to look in more places, so I store the hints from
    #       ENV{TBB_ARCH_PLATFORM} and the TBB_ARCHITECTURE and TBB_COMPILER
    #       variables and search them both.
    set (_TBB_LIBRARY_DIR "${_TBB_INSTALL_DIR}/${_TBB_ARCHITECTURE}/${_TBB_COMPILER}/lib" ${_TBB_LIBRARY_DIR})
endif ((NOT ${TBB_ARCHITECTURE} STREQUAL "") AND (NOT ${TBB_COMPILER} STREQUAL ""))

# GvdB: Mac OS X distribution places libraries directly in lib directory.
list(APPEND _TBB_LIBRARY_DIR ${_TBB_INSTALL_DIR}/lib)

# Jiri: No reason not to check the default paths. From recent versions,
#       tbbvars has started exporting the LIBRARY_PATH and LD_LIBRARY_PATH
#       variables, which now point to the directories of the lib files.
#       It all makes more sense to use the ${_TBB_LIBRARY_DIR} as a HINTS
#       argument instead of the implicit PATHS as it isn't hard-coded
#       but computed by system introspection. Searching the LIBRARY_PATH
#       and LD_LIBRARY_PATH environment variables is now even more important
#       that tbbvars doesn't export TBB_ARCH_PLATFORM and it facilitates
#       the use of TBB built from sources.
#  LR: search first with NO_DEFAULT_PATH...
find_library(TBB_RELEASE_LIBRARY ${_TBB_LIB_RELEASE_NAME} HINTS ${_TBB_LIBRARY_DIR}
        PATHS ENV LIBRARY_PATH ENV LD_LIBRARY_PATH NO_DEFAULT_PATH)
find_library(TBB_MALLOC_RELEASE_LIBRARY ${_TBB_LIB_MALLOC_RELEASE_NAME} HINTS ${_TBB_LIBRARY_DIR}
        PATHS ENV LIBRARY_PATH ENV LD_LIBRARY_PATH NO_DEFAULT_PATH)
find_library(TBB_MALLOCPROXY_RELEASE_LIBRARY ${_TBB_LIB_MALLOCPROXY_RELEASE_NAME} HINTS ${_TBB_LIBRARY_DIR}
        PATHS ENV LIBRARY_PATH ENV LD_LIBRARY_PATH NO_DEFAULT_PATH)
if(NOT TBB_RELEASE_LIBRARY OR NOT TBB_MALLOC_RELEASE_LIBRARY OR NOT TBB_MALLOCPROXY_RELEASE_LIBRARY)
# LR: ... and then search again with NO_DEFAULT_PATH if nothing was found
#  in hinted paths
    find_library(TBB_RELEASE_LIBRARY ${_TBB_LIB_RELEASE_NAME} HINTS ${_TBB_LIBRARY_DIR}
            PATHS ENV LIBRARY_PATH ENV LD_LIBRARY_PATH)
    find_library(TBB_MALLOC_RELEASE_LIBRARY ${_TBB_LIB_MALLOC_RELEASE_NAME} HINTS ${_TBB_LIBRARY_DIR}
            PATHS ENV LIBRARY_PATH ENV LD_LIBRARY_PATH)
    find_library(TBB_MALLOCPROXY_RELEASE_LIBRARY ${_TBB_LIB_MALLOCPROXY_RELEASE_NAME} HINTS ${_TBB_LIBRARY_DIR}
            PATHS ENV LIBRARY_PATH ENV LD_LIBRARY_PATH)
endif()

#Extract path from TBB_RELEASE_LIBRARY name
get_filename_component(TBB_RELEASE_LIBRARY_DIR ${TBB_RELEASE_LIBRARY} PATH)

#TBB_CORRECT_LIB_DIR(TBB_RELEASE_LIBRARY)
#TBB_CORRECT_LIB_DIR(TBB_MALLOC_RELEASE_LIBRARY)
#TBB_CORRECT_LIB_DIR(TBB_MALLOCPROXY_RELEASE_LIBRARY)
mark_as_advanced(TBB_RELEASE_LIBRARY TBB_MALLOC_RELEASE_LIBRARY TBB_MALLOCPROXY_RELEASE_LIBRARY)

#-- Look for debug libraries
# Jiri: Changed the same way as for the release libraries.
find_library(TBB_DEBUG_LIBRARY ${_TBB_LIB_DEBUG_NAME} HINTS ${_TBB_LIBRARY_DIR}
        PATHS ENV LIBRARY_PATH ENV LD_LIBRARY_PATH NO_DEFAULT_PATH)
find_library(TBB_MALLOC_DEBUG_LIBRARY ${_TBB_LIB_MALLOC_DEBUG_NAME} HINTS ${_TBB_LIBRARY_DIR}
        PATHS ENV LIBRARY_PATH ENV LD_LIBRARY_PATH NO_DEFAULT_PATH)
find_library(TBB_MALLOCPROXY_DEBUG_LIBRARY ${_TBB_LIB_MALLOCPROXY_DEBUG_NAME} HINTS ${_TBB_LIBRARY_DIR}
        PATHS ENV LIBRARY_PATH ENV LD_LIBRARY_PATH NO_DEFAULT_PATH)
if(NOT TBB_DEBUG_LIBRARY OR NOT TBB_MALLOC_DEBUG_LIBRARY OR NOT TBB_MALLOCPROXY_DEBUG_LIBRARY)
    find_library(TBB_DEBUG_LIBRARY ${_TBB_LIB_DEBUG_NAME} HINTS ${_TBB_LIBRARY_DIR}
            PATHS ENV LIBRARY_PATH ENV LD_LIBRARY_PATH)
    find_library(TBB_MALLOC_DEBUG_LIBRARY ${_TBB_LIB_MALLOC_DEBUG_NAME} HINTS ${_TBB_LIBRARY_DIR}
            PATHS ENV LIBRARY_PATH ENV LD_LIBRARY_PATH)
    find_library(TBB_MALLOCPROXY_DEBUG_LIBRARY ${_TBB_LIB_MALLOCPROXY_DEBUG_NAME} HINTS ${_TBB_LIBRARY_DIR}
            PATHS ENV LIBRARY_PATH ENV LD_LIBRARY_PATH)
endif()

# Jiri: Self-built TBB stores the debug libraries in a separate directory.
#       Extract path from TBB_DEBUG_LIBRARY name
get_filename_component(TBB_DEBUG_LIBRARY_DIR ${TBB_DEBUG_LIBRARY} PATH)

#TBB_CORRECT_LIB_DIR(TBB_DEBUG_LIBRARY)
#TBB_CORRECT_LIB_DIR(TBB_MALLOC_DEBUG_LIBRARY)
#TBB_CORRECT_LIB_DIR(TBB_MALLOCPROXY_DEBUG_LIBRARY)
mark_as_advanced(TBB_DEBUG_LIBRARY TBB_MALLOC_DEBUG_LIBRARY TBB_MALLOCPROXY_DEBUG_LIBRARY)

if (TBB_INCLUDE_DIR)
    if (TBB_RELEASE_LIBRARY)
        set (TBB_FOUND "YES")

        # NOTE: Removed because we don't want to link with the malloc_proxy by default
        #if (NOT "${TBB_MALLOCPROXY_RELEASE_LIBRARY}" STREQUAL "TBB_MALLOCPROXY_RELEASE_LIBRARY-NOTFOUND")
        #    mark_as_advanced(TBB_MALLOCPROXY_RELEASE_LIBRARY)
        #    set (_TBB_MALLOCPROXY optimized ${TBB_MALLOCPROXY_RELEASE_LIBRARY})
        #endif (NOT "${TBB_MALLOCPROXY_RELEASE_LIBRARY}" STREQUAL "TBB_MALLOCPROXY_RELEASE_LIBRARY-NOTFOUND")
        #if (NOT "${TBB_MALLOCPROXY_DEBUG_LIBRARY}" STREQUAL "TBB_MALLOCPROXY_DEBUG_LIBRARY-NOTFOUND")
        #    mark_as_advanced(TBB_MALLOCPROXY_DEBUG_LIBRARY)
        #    set (_TBB_MALLOCPROXY ${_TBB_MALLOCPROXY} debug ${TBB_MALLOCPROXY_DEBUG_LIBRARY})
        #endif (NOT "${TBB_MALLOCPROXY_DEBUG_LIBRARY}" STREQUAL "TBB_MALLOCPROXY_DEBUG_LIBRARY-NOTFOUND")

        # TBB release library
        set (ALL_TBB_LIBRARIES optimized ${TBB_RELEASE_LIBRARY})

        # TBB debug library found?
        if (TBB_DEBUG_LIBRARY)
            list(APPEND ALL_TBB_LIBRARIES debug ${TBB_DEBUG_LIBRARY})
        else (TBB_DEBUG_LIBRARY)
            # Otherwise, link with the release library even in debug mode
            list(APPEND ALL_TBB_LIBRARIES debug ${TBB_RELEASE_LIBRARY})
        endif (TBB_DEBUG_LIBRARY)

        # TBB malloc - release
        if (TBB_MALLOC_RELEASE_LIBRARY)
            list(APPEND ALL_TBB_LIBRARIES optimized ${TBB_MALLOC_RELEASE_LIBRARY})

            # TBB malloc - debug
            if (TBB_MALLOC_DEBUG_LIBRARY)
                list(APPEND ALL_TBB_LIBRARIES debug ${TBB_MALLOC_DEBUG_LIBRARY})
            else (TBB_MALLOC_DEBUG_LIBRARY)
                list(APPEND ALL_TBB_LIBRARIES debug ${TBB_MALLOC_RELEASE_LIBRARY})
            endif (TBB_MALLOC_DEBUG_LIBRARY)
        endif (TBB_MALLOC_RELEASE_LIBRARY)

        if(UNIX AND NOT APPLE)
            # On Fedora, code using TBB might need -pthread

            # First check without pthread
            try_TBB_with_pthread(TBB_without_pthread)

            if(NOT TBB_without_pthread)
                # Then check with -pthread
                try_TBB_with_pthread(TBB_with_pthread -pthread)
                if(TBB_with_pthread)
                    list(APPEND ALL_TBB_LIBRARIES general -pthread)
                endif(TBB_with_pthread)
                endif(NOT TBB_without_pthread)
        endif(UNIX AND NOT APPLE)

        set (TBB_LIBRARIES ${ALL_TBB_LIBRARIES}
             CACHE PATH "TBB libraries" FORCE)

        # Include dirs
        set (TBB_INCLUDE_DIRS ${TBB_INCLUDE_DIR} CACHE PATH "TBB include directory" FORCE)

        # Library dirs
        if( "${TBB_DEBUG_LIBRARY_DIR}" STREQUAL "" OR "${TBB_RELEASE_LIBRARY_DIR}" STREQUAL "${TBB_DEBUG_LIBRARY_DIR}" )
            set (TBB_LIBRARY_DIRS
                ${TBB_RELEASE_LIBRARY_DIR}
                CACHE PATH "TBB library directories" FORCE)
        else( "${TBB_DEBUG_LIBRARY_DIR}" STREQUAL "" OR "${TBB_RELEASE_LIBRARY_DIR}" STREQUAL "${TBB_DEBUG_LIBRARY_DIR}" )
            set (TBB_LIBRARY_DIRS
                ${TBB_RELEASE_LIBRARY_DIR} ${TBB_DEBUG_LIBRARY_DIR}
                CACHE PATH "TBB library directories" FORCE)
        endif( "${TBB_DEBUG_LIBRARY_DIR}" STREQUAL "" OR "${TBB_RELEASE_LIBRARY_DIR}" STREQUAL "${TBB_DEBUG_LIBRARY_DIR}" )

        message(STATUS "Found Intel TBB")
    endif (TBB_RELEASE_LIBRARY)
endif (TBB_INCLUDE_DIR)

if (NOT TBB_FOUND)
    if(NOT TBB_FIND_QUIETLY)
        message("ERROR: Intel TBB NOT found! Please define the TBBROOT (or TBB_INSTALL_DIR) and/or TBB_ARCH_PLATFORM environment variables.")
        message(STATUS "Looked for Threading Building Blocks in ${_TBB_INSTALL_DIR}")
    endif(NOT TBB_FIND_QUIETLY)
    SET(TBB_INSTALL_DIR "TBB_INSTALL_DIR_NOT_FOUND" CACHE STRING "Intel TBB install directory")
    # do only throw fatal, if this pkg is REQUIRED
    if (TBB_FIND_REQUIRED)
        message(FATAL_ERROR "Could NOT find TBB library.")
    endif (TBB_FIND_REQUIRED)
endif (NOT TBB_FOUND)

endif (NOT _TBB_INSTALL_DIR)

if (TBB_FOUND)
    set(TBB_INTERFACE_VERSION 0)
    FILE(READ "${TBB_INCLUDE_DIRS}/tbb/tbb_stddef.h" _TBB_VERSION_CONTENTS)
    STRING(REGEX REPLACE ".*#define TBB_INTERFACE_VERSION ([0-9]+).*" "\\1" TBB_INTERFACE_VERSION "${_TBB_VERSION_CONTENTS}")
    set(TBB_INTERFACE_VERSION "${TBB_INTERFACE_VERSION}")
endif (TBB_FOUND)

set(TBB_USE_FILE "UseTBB")

### ** Emacs settings **
### Local Variables:
### cmake-tab-width: 4
### End:
