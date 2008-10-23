# Find LAPACK library
#
# This module finds an installed library that implements the LAPACK
# linear-algebra interface (see http://www.netlib.org/lapack/).
# The approach follows mostly that taken for the autoconf macro file, acx_lapack.m4
# (distributed at http://ac-archive.sourceforge.net/ac-archive/acx_lapack.html).
#
# This module sets the following variables:
#  LAPACK_FOUND - set to true if a library implementing the LAPACK interface
#    is found
#  LAPACK_DEFINITIONS - Compilation options to use LAPACK
#  LAPACK_LINKER_FLAGS - Linker flags to use LAPACK (excluding -l
#    and -L).
#  LAPACK_LIBRARIES_DIR - Directories containing the LAPACK libraries. 
#     May be null if LAPACK_LIBRARIES contains libraries name using full path.
#  LAPACK_LIBRARIES - List of libraries to link against LAPACK interface.
#     May be null if the compiler supports auto-link (e.g. VC++).
#
# This module was modified by CGAL team:
# - find LAPACK library shipped with TAUCS
# - find libraries for a C++ compiler, instead of Fortran
# - added LAPACK_DEFINITIONS and LAPACK_LIBRARIES_DIR
# - removed LAPACK95_LIBRARIES
#
# TODO (CGAL): 
# - use a C++ compiler instead of a Fortran one
# - try to be compatible with CMake 2.4
# - find CLAPACK (http://www.netlib.org/clapack)?


# CheckFortranFunctionExists is new in CMake 2.6
CMAKE_MINIMUM_REQUIRED(VERSION 2.6 FATAL_ERROR)

include(CheckFortranFunctionExists)

include(GeneratorSpecificSettings)


# This macro checks for the existence of the combination of fortran libraries
# given by _list.  If the combination is found, this macro checks (using the
# Check_Fortran_Function_Exists macro) whether can link against that library
# combination using the name of a routine given by _name using the linker
# flags given by _flags.  If the combination of libraries is found and passes
# the link test, LIBRARIES is set to the list of complete library paths that
# have been found.  Otherwise, LIBRARIES is set to FALSE.

# N.B. _prefix is the prefix applied to the names of all cached variables that
# are generated internally and marked advanced by this macro.
macro(check_lapack_libraries LIBRARIES _prefix _name _flags _list _blas _path)
  #message("DEBUG: check_lapack_libraries(${_list} in ${_path})")
  set(_libraries_work TRUE)
  set(${LIBRARIES})
  set(_combined_name)
  foreach(_library ${_list})
    set(_combined_name ${_combined_name}_${_library})

    if(_libraries_work)
      if ( WIN32 )
        find_library(${_prefix}_${_library}_LIBRARY
                    NAMES ${_library}
                    PATHS ${_path} ENV LIB
        )
      endif ( WIN32 )
      if ( APPLE )
        find_library(${_prefix}_${_library}_LIBRARY
                    NAMES ${_library}
                    PATHS ${_path} /usr/local/lib /usr/lib /usr/local/lib64 /usr/lib64 ENV DYLD_LIBRARY_PATH
        )
      else ( APPLE )
        find_library(${_prefix}_${_library}_LIBRARY
                    NAMES ${_library}
                    PATHS ${_path} /usr/local/lib /usr/lib /usr/local/lib64 /usr/lib64 ENV LD_LIBRARY_PATH
        )
      endif( APPLE )
      mark_as_advanced(${_prefix}_${_library}_LIBRARY)
      set(${LIBRARIES} ${${LIBRARIES}} ${${_prefix}_${_library}_LIBRARY})
      set(_libraries_work ${${_prefix}_${_library}_LIBRARY})
    endif(_libraries_work)
  endforeach(_library ${_list})

  if(_libraries_work)
    # Test this combination of libraries.
    set(CMAKE_REQUIRED_LIBRARIES ${_flags} ${${LIBRARIES}} ${_blas})
    #message("DEBUG: CMAKE_REQUIRED_LIBRARIES = ${CMAKE_REQUIRED_LIBRARIES}")
    check_fortran_function_exists(${_name} ${_prefix}${_combined_name}_WORKS)
    set(CMAKE_REQUIRED_LIBRARIES)
    mark_as_advanced(${_prefix}${_combined_name}_WORKS)
    set(_libraries_work ${${_prefix}${_combined_name}_WORKS})
  endif(_libraries_work)

  if(NOT _libraries_work)
    set(${LIBRARIES} FALSE)
  endif(NOT _libraries_work)
  #message("DEBUG: ${LIBRARIES} = ${${LIBRARIES}}")
endmacro(check_lapack_libraries)


#
# main
#

# Is it already configured?
if (LAPACK_LIBRARIES_DIR OR LAPACK_LIBRARIES) 

  set(LAPACK_FOUND TRUE)

else(LAPACK_LIBRARIES_DIR OR LAPACK_LIBRARIES)

  # unused (yet)
  set(LAPACK_LINKER_FLAGS)

  # Look first for the LAPACK distributed with CGAL in auxiliary/taucs.
  # Set CGAL_TAUCS_FOUND, CGAL_TAUCS_INCLUDE_DIR and CGAL_TAUCS_LIBRARIES_DIR.
  include(CGAL_Locate_CGAL_TAUCS)

  # Search for LAPACK libraries in ${CGAL_TAUCS_LIBRARIES_DIR} (LAPACK shipped with CGAL),
  # else in $LAPACK_LIB_DIR environment variable.
  if(CGAL_TAUCS_FOUND AND AUTO_LINK_ENABLED)

    # if VC++: done
    #message("DEBUG: LAPACK: VC++ case")
    set( LAPACK_LIBRARIES_DIR  "${CGAL_TAUCS_LIBRARIES_DIR}"
                               CACHE FILEPATH "Directories containing the LAPACK libraries")

  else(CGAL_TAUCS_FOUND AND AUTO_LINK_ENABLED)

    #message("DEBUG: LAPACK: Unix case")

    # If Unix, we need a Fortran compiler
    if (NOT MSVC) # safety: enable_language(Fortran) is broken for VC++
      enable_language(Fortran OPTIONAL)
    endif()

    # LAPACK requires BLAS
    find_package(BLAS QUIET)

    if(CMAKE_Fortran_COMPILER_WORKS AND BLAS_FOUND)

      #intel mkl lapack?
      if(NOT LAPACK_LIBRARIES)
        check_lapack_libraries(
        LAPACK_LIBRARIES
        LAPACK
        cheev
        ""
        "mkl_lapack"
        "${BLAS_LIBRARIES}"
        "${CGAL_TAUCS_LIBRARIES_DIR} $ENV{LAPACK_LIB_DIR}"
        )
        if(LAPACK_LIBRARIES)
          # Use f2c calling convention
          set( LAPACK_DEFINITIONS "-DCGAL_USE_F2C" )
        endif()
      endif()

      #acml lapack?
      if(NOT LAPACK_LIBRARIES)
        check_lapack_libraries(
        LAPACK_LIBRARIES
        LAPACK
        cheev
        ""
        "acml"
        "${BLAS_LIBRARIES}"
        "${CGAL_TAUCS_LIBRARIES_DIR} $ENV{LAPACK_LIB_DIR}"
        )
      endif()

      # Apple LAPACK library?
      if(NOT LAPACK_LIBRARIES)
        check_lapack_libraries(
        LAPACK_LIBRARIES
        LAPACK
        cheev
        ""
        "Accelerate"
        "${BLAS_LIBRARIES}"
        "${CGAL_TAUCS_LIBRARIES_DIR} $ENV{LAPACK_LIB_DIR}"
        )
        if(LAPACK_LIBRARIES)
          # Use f2c calling convention
          set( LAPACK_DEFINITIONS "-DCGAL_USE_F2C" )
        endif()
      endif()

      if ( NOT LAPACK_LIBRARIES )
        check_lapack_libraries(
        LAPACK_LIBRARIES
        LAPACK
        cheev
        ""
        "vecLib"
        "${BLAS_LIBRARIES}"
        "${CGAL_TAUCS_LIBRARIES_DIR} $ENV{LAPACK_LIB_DIR}"
        )
        if(LAPACK_LIBRARIES)
          # Use f2c calling convention
          set( LAPACK_DEFINITIONS "-DCGAL_USE_F2C" )
        endif()
      endif ( NOT LAPACK_LIBRARIES )

      # Generic LAPACK library?
      # This configuration *must* be the last try as this library is notably slow.
      if ( NOT LAPACK_LIBRARIES )
        check_lapack_libraries(
        LAPACK_LIBRARIES
        LAPACK
        cheev
        ""
        "lapack"
        "${BLAS_LIBRARIES}"
        "${CGAL_TAUCS_LIBRARIES_DIR} $ENV{LAPACK_LIB_DIR}"
        )
        if(LAPACK_LIBRARIES)
          # Use f2c calling convention
          set( LAPACK_DEFINITIONS "-DCGAL_USE_F2C" )
        endif()
      endif()

    else(CMAKE_Fortran_COMPILER_WORKS AND BLAS_FOUND)

      message(STATUS "FindLAPACK.cmake requires a Fortran compiler")

    endif(CMAKE_Fortran_COMPILER_WORKS AND BLAS_FOUND)

    # Add variables to cache
    set( LAPACK_DEFINITIONS   "${LAPACK_DEFINITIONS}" 
                              CACHE FILEPATH "Compilation options to use LAPACK" )
    set( LAPACK_LINKER_FLAGS  "${LAPACK_LINKER_FLAGS}" 
                              CACHE FILEPATH "Linker flags to use LAPACK" )
    set( LAPACK_LIBRARIES     "${LAPACK_LIBRARIES}" 
                              CACHE FILEPATH "LAPACK libraries name" )

  endif(CGAL_TAUCS_FOUND AND AUTO_LINK_ENABLED)

  if(LAPACK_LIBRARIES_DIR OR LAPACK_LIBRARIES)
    set(LAPACK_FOUND TRUE)
  else()
    set(LAPACK_FOUND FALSE)
  endif()

  if(NOT LAPACK_FIND_QUIETLY)
    if(LAPACK_FOUND)
      message(STATUS "A library with LAPACK API found.")
    else(LAPACK_FOUND)
      if(LAPACK_FIND_REQUIRED)
        message(FATAL_ERROR "A required library with LAPACK API not found. Please specify library location.")
      else()
        message(STATUS "A library with LAPACK API not found. Please specify library location.")
      endif()
    endif(LAPACK_FOUND)
  endif(NOT LAPACK_FIND_QUIETLY)

  #message("DEBUG: LAPACK_DEFINITIONS = ${LAPACK_DEFINITIONS}")
  #message("DEBUG: LAPACK_LINKER_FLAGS = ${LAPACK_LINKER_FLAGS}")
  #message("DEBUG: LAPACK_LIBRARIES = ${LAPACK_LIBRARIES}")
  #message("DEBUG: LAPACK_LIBRARIES_DIR = ${LAPACK_LIBRARIES_DIR}")
  #message("DEBUG: LAPACK_FOUND = ${LAPACK_FOUND}")

endif(LAPACK_LIBRARIES_DIR OR LAPACK_LIBRARIES)

#mark_as_advanced(LAPACK_DEFINITIONS)
#mark_as_advanced(LAPACK_LINKER_FLAGS)
#mark_as_advanced(LAPACK_LIBRARIES)
#mark_as_advanced(LAPACK_LIBRARIES_DIR)
