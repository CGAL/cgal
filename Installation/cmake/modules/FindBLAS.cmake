# Find BLAS library
#
# This module finds an installed library that implements the BLAS
# linear-algebra interface (see http://www.netlib.org/blas/).
# The list of libraries searched for is mainly taken
# from the autoconf macro file, acx_blas.m4 (distributed at
# http://ac-archive.sourceforge.net/ac-archive/acx_blas.html).
#
# This module sets the following variables:
#  BLAS_FOUND - set to true if a library implementing the BLAS interface
#    is found
#  BLAS_DEFINITIONS - Compilation options to use BLAS
#  BLAS_LINKER_FLAGS - Linker flags to use BLAS (excluding -l
#    and -L).
#  BLAS_LIBRARIES_DIR - Directories containing the BLAS libraries.
#     May be null if BLAS_LIBRARIES contains libraries name using full path.
#  BLAS_LIBRARIES - List of libraries to link against BLAS interface.
#     May be null if the compiler supports auto-link (e.g. VC++).
#
# This module was modified by CGAL team:
# - find BLAS library shipped with TAUCS
# - find libraries for a C++ compiler, instead of Fortran
# - added BLAS_DEFINITIONS and BLAS_LIBRARIES_DIR
# - removed BLAS95_LIBRARIES
#
# TODO (CGAL):
# - use a C++ compiler instead of a Fortran one
# - try to be compatible with CMake 2.4
# - find CBLAS (http://www.netlib.org/cblas)?


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
macro(check_fortran_libraries LIBRARIES _prefix _name _flags _list _path)
  #message("DEBUG: check_fortran_libraries(${_list} in ${_path})")
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
    set(CMAKE_REQUIRED_LIBRARIES ${_flags} ${${LIBRARIES}})
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
endmacro(check_fortran_libraries)


# This macro setup BLAS variables to link with Fortran libraries:
# - it adds "-DCGAL_USE_F2C" to use f2c calling convention
# - its links with the f2c library if needed
macro(append_f2c)
  # Use f2c calling convention
  set( BLAS_DEFINITIONS  ${BLAS_DEFINITIONS} "-DCGAL_USE_F2C" )

  # Some C++ linkers require f2c to link with Fortran libraries.
  # Implementation note: I do not know which ones, thus I just add
  #                      the f2c library if it is available.
  find_package( F2C QUIET )
  if ( F2C_FOUND )
    set( BLAS_DEFINITIONS  ${BLAS_DEFINITIONS} ${F2C_DEFINITIONS} )
    set( BLAS_LIBRARIES    ${BLAS_LIBRARIES} ${F2C_LIBRARIES} )
  endif()
endmacro(append_f2c)


#
# main
#

# Is it already configured?
if (BLAS_LIBRARIES_DIR OR BLAS_LIBRARIES)

  set(BLAS_FOUND TRUE)

else(BLAS_LIBRARIES_DIR OR BLAS_LIBRARIES)

  # Look first for the BLAS distributed with CGAL in auxiliary/taucs.
  # Set CGAL_TAUCS_FOUND, CGAL_TAUCS_INCLUDE_DIR and CGAL_TAUCS_LIBRARIES_DIR.
  include(CGAL_Locate_CGAL_TAUCS)

  # Search for BLAS libraries in ${CGAL_TAUCS_LIBRARIES_DIR} (BLAS shipped with CGAL),
  # else in $BLAS_LIB_DIR environment variable.
  if(CGAL_TAUCS_FOUND AND CGAL_AUTO_LINK_ENABLED)

    # if VC++: done
    #message("DEBUG: BLAS: VC++ case")
    set( BLAS_LIBRARIES_DIR  "${CGAL_TAUCS_LIBRARIES_DIR}"
                             CACHE FILEPATH "Directories containing the BLAS libraries")

  else(CGAL_TAUCS_FOUND AND CGAL_AUTO_LINK_ENABLED)

    #message("DEBUG: BLAS: Unix case")

    # If Unix, we need a Fortran compiler
    if (NOT MSVC) # safety: enable_language(Fortran) is broken for VC++
      enable_language(Fortran OPTIONAL)
    endif()
    if (CMAKE_Fortran_COMPILER_WORKS)

      # BLAS in ATLAS library? (http://math-atlas.sourceforge.net/)
      if(NOT BLAS_LIBRARIES)
        check_fortran_libraries(
        BLAS_LIBRARIES
        BLAS
        cblas_dgemm
        ""
        "cblas;f77blas;atlas"
        "${CGAL_TAUCS_LIBRARIES_DIR} $ENV{BLAS_LIB_DIR}"
        )
        if(BLAS_LIBRARIES)
          # Some C++ linkers require f2c to link with Fortran libraries
          append_f2c()
        endif()
      endif()

      # BLAS in PhiPACK libraries? (requires generic BLAS lib, too)
      if(NOT BLAS_LIBRARIES)
        check_fortran_libraries(
        BLAS_LIBRARIES
        BLAS
        sgemm
        ""
        "sgemm;dgemm;blas"
        "${CGAL_TAUCS_LIBRARIES_DIR} $ENV{BLAS_LIB_DIR}"
        )
        if(BLAS_LIBRARIES)
          # Some C++ linkers require f2c to link with Fortran libraries
          append_f2c()
        endif()
      endif()

      # BLAS in Alpha CXML library?
      if(NOT BLAS_LIBRARIES)
        check_fortran_libraries(
        BLAS_LIBRARIES
        BLAS
        sgemm
        ""
        "cxml"
        "${CGAL_TAUCS_LIBRARIES_DIR} $ENV{BLAS_LIB_DIR}"
        )
      endif()

      # BLAS in Alpha DXML library? (now called CXML, see above)
      if(NOT BLAS_LIBRARIES)
        check_fortran_libraries(
        BLAS_LIBRARIES
        BLAS
        sgemm
        ""
        "dxml"
        "${CGAL_TAUCS_LIBRARIES_DIR} $ENV{BLAS_LIB_DIR}"
        )
      endif()

      # BLAS in Sun Performance library?
      if(NOT BLAS_LIBRARIES)
        check_fortran_libraries(
        BLAS_LIBRARIES
        BLAS
        sgemm
        "-xlic_lib=sunperf"
        "sunperf;sunmath"
        "${CGAL_TAUCS_LIBRARIES_DIR} $ENV{BLAS_LIB_DIR}"
        )
        if(BLAS_LIBRARIES)
          # Use f2c calling convention
          set( BLAS_DEFINITIONS  "-DCGAL_USE_F2C" )
          # Extra library
          set(BLAS_LINKER_FLAGS "-xlic_lib=sunperf")
        endif()
      endif()

      # BLAS in SCSL library?  (SGI/Cray Scientific Library)
      if(NOT BLAS_LIBRARIES)
        check_fortran_libraries(
        BLAS_LIBRARIES
        BLAS
        sgemm
        ""
        "scsl"
        "${CGAL_TAUCS_LIBRARIES_DIR} $ENV{BLAS_LIB_DIR}"
        )
      endif()

      # BLAS in SGIMATH library?
      if(NOT BLAS_LIBRARIES)
        check_fortran_libraries(
        BLAS_LIBRARIES
        BLAS
        sgemm
        ""
        "complib.sgimath"
        "${CGAL_TAUCS_LIBRARIES_DIR} $ENV{BLAS_LIB_DIR}"
        )
      endif()

      # BLAS in IBM ESSL library? (requires generic BLAS lib, too)
      if(NOT BLAS_LIBRARIES)
        check_fortran_libraries(
        BLAS_LIBRARIES
        BLAS
        sgemm
        ""
        "essl;blas"
        "${CGAL_TAUCS_LIBRARIES_DIR} $ENV{BLAS_LIB_DIR}"
        )
      endif()

      #BLAS in intel mkl 10 library? (em64t 64bit)
      if(NOT BLAS_LIBRARIES)
        check_fortran_libraries(
        BLAS_LIBRARIES
        BLAS
        sgemm
        ""
        "mkl_intel_lp64;mkl_intel_thread;mkl_core;guide;pthread"
        "${CGAL_TAUCS_LIBRARIES_DIR} $ENV{BLAS_LIB_DIR}"
        )
        if(BLAS_LIBRARIES)
          # Use f2c calling convention
          set( BLAS_DEFINITIONS  "-DCGAL_USE_F2C" )
        endif()
      endif()

      ### windows version of intel mkl 10?
      if(NOT BLAS_LIBRARIES)
        check_fortran_libraries(
        BLAS_LIBRARIES
        BLAS
        SGEMM
        ""
        "mkl_c_dll;mkl_intel_thread_dll;mkl_core_dll;libguide40"
        "${CGAL_TAUCS_LIBRARIES_DIR} $ENV{BLAS_LIB_DIR}"
        )
        if(BLAS_LIBRARIES)
          # Use f2c calling convention
          set( BLAS_DEFINITIONS  "-DCGAL_USE_F2C" )
        endif()
      endif()

      #older versions of intel mkl libs

      # BLAS in intel mkl library? (shared)
      if(NOT BLAS_LIBRARIES)
        check_fortran_libraries(
        BLAS_LIBRARIES
        BLAS
        sgemm
        ""
        "mkl;guide;pthread"
        "${CGAL_TAUCS_LIBRARIES_DIR} $ENV{BLAS_LIB_DIR}"
        )
        if(BLAS_LIBRARIES)
          # Use f2c calling convention
          set( BLAS_DEFINITIONS  "-DCGAL_USE_F2C" )
        endif()
      endif()

      #BLAS in intel mkl library? (static, 32bit)
      if(NOT BLAS_LIBRARIES)
        check_fortran_libraries(
        BLAS_LIBRARIES
        BLAS
        sgemm
        ""
        "mkl_ia32;guide;pthread"
        "${CGAL_TAUCS_LIBRARIES_DIR} $ENV{BLAS_LIB_DIR}"
        )
        if(BLAS_LIBRARIES)
          # Use f2c calling convention
          set( BLAS_DEFINITIONS  "-DCGAL_USE_F2C" )
        endif()
      endif()

      #BLAS in intel mkl library? (static, em64t 64bit)
      if(NOT BLAS_LIBRARIES)
        check_fortran_libraries(
        BLAS_LIBRARIES
        BLAS
        sgemm
        ""
        "mkl_em64t;guide;pthread"
        "${CGAL_TAUCS_LIBRARIES_DIR} $ENV{BLAS_LIB_DIR}"
        )
        if(BLAS_LIBRARIES)
          # Use f2c calling convention
          set( BLAS_DEFINITIONS  "-DCGAL_USE_F2C" )
        endif()
      endif()

      #BLAS in acml library?
      if(NOT BLAS_LIBRARIES)
        check_fortran_libraries(
        BLAS_LIBRARIES
        BLAS
        sgemm
        ""
        "acml"
        "${CGAL_TAUCS_LIBRARIES_DIR} $ENV{BLAS_LIB_DIR}"
        )
      endif()

      # Apple BLAS library?
      if(NOT BLAS_LIBRARIES)
        check_fortran_libraries(
        BLAS_LIBRARIES
        BLAS
        cblas_dgemm
        ""
        "Accelerate"
        "${CGAL_TAUCS_LIBRARIES_DIR} $ENV{BLAS_LIB_DIR}"
        )
        if(BLAS_LIBRARIES)
          # Use f2c calling convention
          set( BLAS_DEFINITIONS  "-DCGAL_USE_F2C" )
        endif()
      endif()

      if ( NOT BLAS_LIBRARIES )
        check_fortran_libraries(
        BLAS_LIBRARIES
        BLAS
        cblas_dgemm
        ""
        "vecLib"
        "${CGAL_TAUCS_LIBRARIES_DIR} $ENV{BLAS_LIB_DIR}"
        )
        if(BLAS_LIBRARIES)
          # Use f2c calling convention
          set( BLAS_DEFINITIONS  "-DCGAL_USE_F2C" )
        endif()
      endif ( NOT BLAS_LIBRARIES )

      # Generic BLAS library?
      # This configuration *must* be the last try as this library is notably slow.
      if(NOT BLAS_LIBRARIES)
        check_fortran_libraries(
        BLAS_LIBRARIES
        BLAS
        sgemm
        ""
        "blas"
        "${CGAL_TAUCS_LIBRARIES_DIR} $ENV{BLAS_LIB_DIR}"
        )
        if(BLAS_LIBRARIES)
          # Some C++ linkers require f2c to link with Fortran libraries
          append_f2c()
        endif()
      endif()

    else(CMAKE_Fortran_COMPILER_WORKS)

      message(STATUS "FindBLAS.cmake requires a Fortran compiler")

    endif(CMAKE_Fortran_COMPILER_WORKS)

    # Add variables to cache
    set( BLAS_DEFINITIONS   "${BLAS_DEFINITIONS}"
                            CACHE FILEPATH "Compilation options to use BLAS" )
    set( BLAS_LINKER_FLAGS  "${BLAS_LINKER_FLAGS}"
                            CACHE FILEPATH "Linker flags to use BLAS" )
    set( BLAS_LIBRARIES     "${BLAS_LIBRARIES}"
                            CACHE FILEPATH "BLAS libraries name" )

  endif(CGAL_TAUCS_FOUND AND CGAL_AUTO_LINK_ENABLED)

  if(BLAS_LIBRARIES_DIR OR BLAS_LIBRARIES)
    set(BLAS_FOUND TRUE)
  else()
    set(BLAS_FOUND FALSE)
  endif()

  if(NOT BLAS_FIND_QUIETLY)
    if(BLAS_FOUND)
      message(STATUS "A library with BLAS API found.")
    else(BLAS_FOUND)
      if(BLAS_FIND_REQUIRED)
        message(FATAL_ERROR "A required library with BLAS API not found. Please specify library location.")
      else()
        message(STATUS "A library with BLAS API not found. Please specify library location.")
      endif()
    endif(BLAS_FOUND)
  endif(NOT BLAS_FIND_QUIETLY)

  #message("DEBUG: BLAS_DEFINITIONS = ${BLAS_DEFINITIONS}")
  #message("DEBUG: BLAS_LINKER_FLAGS = ${BLAS_LINKER_FLAGS}")
  #message("DEBUG: BLAS_LIBRARIES = ${BLAS_LIBRARIES}")
  #message("DEBUG: BLAS_LIBRARIES_DIR = ${BLAS_LIBRARIES_DIR}")
  #message("DEBUG: BLAS_FOUND = ${BLAS_FOUND}")

endif(BLAS_LIBRARIES_DIR OR BLAS_LIBRARIES)

#mark_as_advanced(BLAS_DEFINITIONS)
#mark_as_advanced(BLAS_LINKER_FLAGS)
#mark_as_advanced(BLAS_LIBRARIES)
#mark_as_advanced(BLAS_LIBRARIES_DIR)
