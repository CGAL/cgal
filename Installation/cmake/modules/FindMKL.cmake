# Find MKL library
#
# This module finds Intel Math Kernel Library (including PARDISO).
# It is based on FindBLAS.cmake.
#
# This module sets the following variables:
#  MKL_FOUND - set to true if MKL is found
#  MKL_INCLUDE_DIR - Directories containing the MKL header files
#  MKL_DEFINITIONS - Compilation options to use MKL
#  MKL_LIBRARIES - List of libraries to link against MKL interface.
#     May be null if the compiler supports auto-link (e.g. VC++).
#  MKL_USE_FILE - The name of the cmake module to include to compile
#     applications or libraries using MKL.


include(CheckFunctionExists)

include(${CMAKE_CURRENT_LIST_DIR}/CGAL_GeneratorSpecificSettings.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/CGAL_Macros.cmake)


# This macro checks for the existence of the combination of fortran libraries
# given by _list.  If the combination is found, this macro checks (using the
# check_function_exists macro) whether can link against that library
# combination using the name of a routine given by _name using the linker
# flags given by _flags.  If the combination of libraries is found and passes
# the link test, LIBRARIES is set to the list of complete library paths that
# have been found and DEFINITIONS to the required definitions.
# Otherwise, LIBRARIES is set to FALSE.
# N.B. _prefix is the prefix applied to the names of all cached variables that
# are generated internally and marked advanced by this macro.
macro(check_fortran_libraries DEFINITIONS LIBRARIES _prefix _name _flags _list _path1 _path2)
  #message("DEBUG: check_fortran_libraries(${_list} in ${_path1} ${_path2} )")

  # Check for the existence of the libraries given by _list
  set(_libraries_found TRUE)
  set(_libraries_work FALSE)
  set(${DEFINITIONS} "")
  set(${LIBRARIES} "")
  set(_combined_name)
  foreach(_library ${_list})
    set(_combined_name ${_combined_name}_${_library})

    if(_libraries_found)
      # search first in ${_path1} ${_path2} 
      find_library(${_prefix}_${_library}_LIBRARY
                  NAMES ${_library}
                  PATHS ${_path1} ${_path2}  NO_DEFAULT_PATH
                  )
      # if not found, search in environment variables and system
      if ( WIN32 )
        find_library(${_prefix}_${_library}_LIBRARY
                    NAMES ${_library}
                    PATHS ENV LIB
                    )
      elseif ( APPLE )
        find_library(${_prefix}_${_library}_LIBRARY
                    NAMES ${_library}
                    PATHS /usr/local/lib /usr/lib /usr/local/lib64 /usr/lib64 ENV DYLD_LIBRARY_PATH
                    )
      else ()
        find_library(${_prefix}_${_library}_LIBRARY
                    NAMES ${_library}
                    PATHS /usr/local/lib /usr/lib /usr/local/lib64 /usr/lib64 ENV LD_LIBRARY_PATH
                    )
      endif()
      #message("DEBUG: find_library(${_library}) = ${${_prefix}_${_library}_LIBRARY}")
      mark_as_advanced(${_prefix}_${_library}_LIBRARY)
      set(${LIBRARIES} ${${LIBRARIES}} ${${_prefix}_${_library}_LIBRARY})
      set(_libraries_found ${${_prefix}_${_library}_LIBRARY})
    endif(_libraries_found)
  endforeach(_library ${_list})
  if(_libraries_found)
    set(_libraries_found ${${LIBRARIES}})
  endif()

  # Test this combination of libraries with the Fortran/f2c interface.
  # We test the Fortran interface first as it is well standardized.
  if(_libraries_found AND NOT _libraries_work)
    set(${DEFINITIONS}  "-D${_prefix}_USE_F2C")
    set(${LIBRARIES}    ${_libraries_found})
    # Some C++ linkers require the f2c library to link with Fortran libraries.
    # I do not know which ones, thus I just add the f2c library if it is available.
    find_package( F2C QUIET )
    if ( F2C_FOUND )
      set(${DEFINITIONS}  ${${DEFINITIONS}} ${F2C_DEFINITIONS})
      set(${LIBRARIES}    ${${LIBRARIES}} ${F2C_LIBRARIES})
    endif()
    set(CMAKE_REQUIRED_DEFINITIONS  ${${DEFINITIONS}})
    set(CMAKE_REQUIRED_LIBRARIES    ${_flags} ${${LIBRARIES}})
    #message("DEBUG: CMAKE_REQUIRED_DEFINITIONS = ${CMAKE_REQUIRED_DEFINITIONS}")
    #message("DEBUG: CMAKE_REQUIRED_LIBRARIES = ${CMAKE_REQUIRED_LIBRARIES}")
    # Check if function exists with f2c calling convention (ie a trailing underscore)
    check_function_exists(${_name}_ ${_prefix}_${_name}_${_combined_name}_f2c_WORKS)
    #message("DEBUG: check_function_exists(${_name}_) = ${${_prefix}_${_name}_${_combined_name}_f2c_WORKS}")
    set(CMAKE_REQUIRED_DEFINITIONS} "")
    set(CMAKE_REQUIRED_LIBRARIES    "")
    mark_as_advanced(${_prefix}_${_name}_${_combined_name}_f2c_WORKS)
    set(_libraries_work ${${_prefix}_${_name}_${_combined_name}_f2c_WORKS})
  endif(_libraries_found AND NOT _libraries_work)

  # If not found, test this combination of libraries with a C interface.
  # A few implementations (ie ACML) provide a C interface. Unfortunately, there is no standard.
  if(_libraries_found AND NOT _libraries_work)
    set(${DEFINITIONS} "")
    set(${LIBRARIES}   ${_libraries_found})
    set(CMAKE_REQUIRED_DEFINITIONS "")
    set(CMAKE_REQUIRED_LIBRARIES   ${_flags} ${${LIBRARIES}})
    #message("DEBUG: CMAKE_REQUIRED_LIBRARIES = ${CMAKE_REQUIRED_LIBRARIES}")
    check_function_exists(${_name} ${_prefix}_${_name}${_combined_name}_WORKS)
    #message("DEBUG: check_function_exists(${_name}) = ${${_prefix}_${_name}${_combined_name}_WORKS}")
    set(CMAKE_REQUIRED_LIBRARIES "")
    mark_as_advanced(${_prefix}_${_name}${_combined_name}_WORKS)
    set(_libraries_work ${${_prefix}_${_name}${_combined_name}_WORKS})
  endif(_libraries_found AND NOT _libraries_work)

  # on failure
  if(NOT _libraries_work)
    set(${DEFINITIONS} "")
    set(${LIBRARIES}   FALSE)
  endif()
  #message("DEBUG: ${DEFINITIONS} = ${${DEFINITIONS}}")
  #message("DEBUG: ${LIBRARIES} = ${${LIBRARIES}}")
endmacro(check_fortran_libraries)


#
# main
#

# Is it already configured?
if (MKL_INCLUDE_DIR AND MKL_LIBRARIES)

  set(MKL_FOUND TRUE)

else()

  ## reset variables
  #set( MKL_INCLUDE_DIR "" )
  #set( MKL_DEFINITIONS "" )
  #set( MKL_LIBRARIES "" )

  # Search for MKL headers
  # in $MKL_INC_DIR environment variable and in usual places.
  find_path(MKL_INCLUDE_DIR
            NAMES mkl.h
            PATHS ENV MKL_INC_DIR NO_DEFAULT_PATH
           )
  find_path(MKL_INCLUDE_DIR
            NAMES mkl.h
           )

  #
  # Search for MKL in possible libraries
  # in $MKL_LIB_DIR environment variable and in usual places.
  #

  # Read environment variables
  fetch_env_var(MKL_LIB_DIR)
  fetch_env_var(INTEL_RTL_LIB_DIR)

  # intel mkl 10 library?
  # TODO: add shared variants
  if (WIN32)
    # intel mkl library? (static, 32bit)
    if(NOT MKL_LIBRARIES)
      check_fortran_libraries(
      MKL_DEFINITIONS
      MKL_LIBRARIES
      MKL
      SGEMM
      ""
      "mkl_solver;mkl_intel_c;mkl_intel_thread;mkl_core;libiomp5md"
      #"mkl_solver_sequential;mkl_intel_c;mkl_sequential;mkl_core"
      "${MKL_LIB_DIR}" "${INTEL_RTL_LIB_DIR}"
      )
    endif()

    # intel mkl library? (static, ia64 and em64t 64 bit)
    if(NOT MKL_LIBRARIES)
      check_fortran_libraries(
      MKL_DEFINITIONS
      MKL_LIBRARIES
      MKL
      SGEMM
      ""
      "mkl_solver_lp64;mkl_intel_lp64;mkl_intel_thread;mkl_core;libiomp5md"
      #"mkl_solver_ilp64_sequential;mkl_intel_ilp64;mkl_sequential;mkl_core"
      "${MKL_LIB_DIR}" "${INTEL_RTL_LIB_DIR}"
      )
    endif()
  else(WIN32)
    # intel mkl library? (static, 32bit)
    if(NOT MKL_LIBRARIES)
      check_fortran_libraries(
      MKL_DEFINITIONS
      MKL_LIBRARIES
      MKL
      sgemm
      ""
      "mkl_solver;mkl_intel;mkl_intel_thread;mkl_core;iomp5;pthread;m"
      "${MKL_LIB_DIR}" "${INTEL_RTL_LIB_DIR}"
      )
    endif()

    # intel mkl library? (static, ia64 and em64t 64 bit)
    if(NOT MKL_LIBRARIES)
      check_fortran_libraries(
      MKL_DEFINITIONS
      MKL_LIBRARIES
      MKL
      sgemm
      ""
      "mkl_solver_lp64;mkl_intel_lp64;mkl_intel_thread;mkl_core;iomp5;pthread;m"
      "${MKL_LIB_DIR}" "${INTEL_RTL_LIB_DIR}"
      )
    endif()
  endif (WIN32)

  # older versions of intel mkl libs

  # intel mkl library? (shared)
  if(NOT MKL_LIBRARIES)
    check_fortran_libraries(
    MKL_DEFINITIONS
    MKL_LIBRARIES
    MKL
    sgemm
    ""
    "mkl_solver;mkl_lapack;mkl;guide;pthread"
    "${MKL_LIB_DIR}" "${INTEL_RTL_LIB_DIR}"
    )
  endif()

  # intel mkl library? (static, 32bit)
  if(NOT MKL_LIBRARIES)
    check_fortran_libraries(
    MKL_DEFINITIONS
    MKL_LIBRARIES
    MKL
    sgemm
    ""
    "mkl_solver;mkl_lapack;mkl_ia32;guide;pthread"
    "${MKL_LIB_DIR}" "${INTEL_RTL_LIB_DIR}"
    )
  endif()

  # intel mkl library? (static, ia64 64bit)
  if(NOT MKL_LIBRARIES)
    check_fortran_libraries(
    MKL_DEFINITIONS
    MKL_LIBRARIES
    MKL
    sgemm
    ""
    "mkl_solver;mkl_lapack;mkl_ipf;guide;pthread"
    "${MKL_LIB_DIR}" "${INTEL_RTL_LIB_DIR}"
    )
  endif()

  # intel mkl library? (static, em64t 64bit)
  if(NOT MKL_LIBRARIES)
    check_fortran_libraries(
    MKL_DEFINITIONS
    MKL_LIBRARIES
    MKL
    sgemm
    ""
    "mkl_solver;mkl_lapack;mkl_em64t;guide;pthread"
    "${MKL_LIB_DIR}" "${INTEL_RTL_LIB_DIR}"
    )
  endif()

  if (MKL_INCLUDE_DIR AND MKL_LIBRARIES)
    set(MKL_FOUND TRUE)
  else()
    set(MKL_FOUND FALSE)
  endif()

  if(NOT MKL_FIND_QUIETLY)
    if(MKL_FOUND)
      message(STATUS "MKL library found.")
    else(MKL_FOUND)
      if(MKL_FIND_REQUIRED)
        message(FATAL_ERROR "MKL is required and not found. Please specify library location.")
      else()
        message(STATUS "MKL not found. Please specify library location.")
      endif()
    endif(MKL_FOUND)
  endif(NOT MKL_FIND_QUIETLY)

  # Add variables to cache
  set( MKL_INCLUDE_DIR   "${MKL_INCLUDE_DIR}"
                          CACHE PATH "Directories containing the MKL header files" FORCE )
  set( MKL_DEFINITIONS   "${MKL_DEFINITIONS}"
                          CACHE STRING "Compilation options to use MKL" FORCE )
  set( MKL_LIBRARIES     "${MKL_LIBRARIES}"
                          CACHE FILEPATH "MKL libraries name" FORCE )

  #message("DEBUG: MKL_INCLUDE_DIR = ${MKL_INCLUDE_DIR}")
  #message("DEBUG: MKL_DEFINITIONS = ${MKL_DEFINITIONS}")
  #message("DEBUG: MKL_LIBRARIES = ${MKL_LIBRARIES}")
  #message("DEBUG: MKL_FOUND = ${MKL_FOUND}")

endif(MKL_INCLUDE_DIR AND MKL_LIBRARIES)

if(MKL_FOUND)
  set(MKL_USE_FILE "CGAL_UseMKL")
endif(MKL_FOUND)
