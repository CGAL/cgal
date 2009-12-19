# Top level CMakeLists.txt for CGAL
# The name of our project is "CGAL".  CMakeLists files in this project can
# refer to the root source directory of the project as ${CMAKE_SOURCE_DIR} or
# ${CMAKE_SOURCE_DIR} and to the root binary directory of the project as
# ${CMAKE_BINARY_DIR} or ${CMAKE_BINARY_DIR}.
project(CGAL)

CMAKE_MINIMUM_REQUIRED(VERSION 2.4.5 FATAL_ERROR)

# This allows else(), endif(), etc... (without repeating the expression)
set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

if ( COMMAND cmake_policy )
  cmake_policy( SET CMP0003 NEW )  
endif( COMMAND cmake_policy )

#--------------------------------------------------------------------------------------------------
#
#                                    -= BASIC SETUP =-
#
#--------------------------------------------------------------------------------------------------


# 
# Read and parse CGAL version number from VERSION file
#

file( READ "VERSION" CGAL_VERSION )
string( REPLACE "-" "." CGAL_VERSION_TOKENS1 ${CGAL_VERSION} )
string( REPLACE "." ";" CGAL_VERSION_TOKENS ${CGAL_VERSION_TOKENS1} )
list( LENGTH CGAL_VERSION_TOKENS CGAL_VERSION_TOKENS_LEN )
list( GET CGAL_VERSION_TOKENS 0 CGAL_MAJOR_VERSION )

set(CGAL_MINOR_VERSION)
set(CGAL_BUILD_VERSION)

if ( CGAL_VERSION_TOKENS_LEN GREATER 1 )
  list( GET CGAL_VERSION_TOKENS 1 CGAL_MINOR_VERSION )
endif()

if ( CGAL_VERSION_TOKENS_LEN GREATER 3 )
  list( GET CGAL_VERSION_TOKENS 3 CGAL_BUILD_VERSION )
elseif( CGAL_VERSION_TOKENS_LEN GREATER 2 )
  list( GET CGAL_VERSION_TOKENS 2 CGAL_BUILD_VERSION )
endif()

message( STATUS "CGAL_MAJOR_VERSION=${CGAL_MAJOR_VERSION}" )
message( STATUS "CGAL_MINOR_VERSION=${CGAL_MINOR_VERSION}" )
message( STATUS "CGAL_BUILD_VERSION=${CGAL_BUILD_VERSION}" )

set( CGAL_SONAME_VERSION "4"     CACHE STRING "soname-version" )
set( CGAL_SOVERSION      "4.0.1" CACHE STRING "soversion" )

mark_as_advanced(CGAL_SONAME_VERSION)
mark_as_advanced(CGAL_SOVERSION)

message( STATUS "CGAL_SONAME_VERSION=${CGAL_SONAME_VERSION}" )
message( STATUS "CGAL_SOVERSION     =${CGAL_SOVERSION}" )

set( CGAL_BUILDING_LIBS TRUE )

set( CGAL_VERSION_DIR     CGAL-${CGAL_VERSION} )
set( CGAL_MODULES_REL_DIR cmake/modules )
set( CGAL_MODULES_DIR     ${CMAKE_SOURCE_DIR}/${CGAL_MODULES_REL_DIR} )

# Where to look first for cmake modules, before ${CMAKE_ROOT}/Modules/ is checked
set(CGAL_CMAKE_MODULE_PATH ${CGAL_MODULES_DIR} )


set(ORIGINAL_CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} )

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CGAL_CMAKE_MODULE_PATH} )


message( STATUS "CGAL_REFERENCE_CACHE_DIR=${CGAL_REFERENCE_CACHE_DIR}" )

if ( RUNNING_CGAL_AUTO_TEST )
  if ( NOT "${CGAL_REFERENCE_CACHE_DIR}" STREQUAL "" )
    if ( EXISTS ${CGAL_REFERENCE_CACHE_DIR} )
      if ( EXISTS ${CGAL_REFERENCE_CACHE_DIR}/CMakeCache.txt )
        message( STATUS "Loading reference cache from ${CGAL_REFERENCE_CACHE_DIR}" )
        load_cache( ${CGAL_REFERENCE_CACHE_DIR} )
      endif()
    endif()
  endif()
endif()

include(CGAL_Common)
include(CGAL_GeneratorSpecificSettings)
include(CGAL_CheckCXXFileRuns)

mark_as_advanced(CMAKE_BACKWARDS_COMPATIBILITY)

message( STATUS "USING CMake version: ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}.${CMAKE_PATCH_VERSION}" )
message( STATUS "System: ${CMAKE_SYSTEM_NAME}" )

# The base name of all CGAL libraries (including "CGAL" itself. That
# variables can be used by foreach() loops that deal with library variables.
set(ALL_CGAL_LIBRARIES CGAL CGAL_Core CGAL_ImageIO CGAL_PDB CGAL_Qt3 CGAL_Qt4 CACHE INTERNAL "Name of all CGAL libraries.")

#--------------------------------------------------------------------------------------------------
#
#                                    -= FLAGS =-
#
#--------------------------------------------------------------------------------------------------

if( MSVC )
  
  uniquely_add_flags ( CGAL_CXX_FLAGS "-D_CRT_SECURE_NO_DEPRECATE;-D_SCL_SECURE_NO_DEPRECATE;-D_CRT_SECURE_NO_WARNINGS;-D_SCL_SECURE_NO_WARNINGS" )
  uniquely_add_flags ( CGAL_CXX_FLAGS "/fp:strict" )
  uniquely_add_flags ( CGAL_CXX_FLAGS "/fp:except-" )

  if ( RUNNING_CGAL_AUTO_TEST )
    set(CMAKE_CXX_WARNING_LEVEL 2 CACHE STRING "MSVC C++ compiler warning level" FORCE)
    mark_as_advanced(CMAKE_CXX_WARNING_LEVEL)
  endif()
    
endif()

if( "${CMAKE_CXX_COMPILER_ID}" MATCHES SunPro )
  message( STATUS "Using SunPro compiler, using STLPort 4." )
  
  uniquely_add_flags( CGAL_CXX_FLAGS "-features=extensions;-library=stlport4;-D_GNU_SOURCE" )
  uniquely_add_flags( CGAL_SHARED_LINKER_FLAGS "-library=stlport4" )
  uniquely_add_flags( CGAL_EXE_LINKER_FLAGS    "-library=stlport4" )
endif()

if( "${CMAKE_CXX_COMPILER}" MATCHES "icl" OR "${CMAKE_CXX_COMPILER}" MATCHES "icpc")
  message( STATUS "Intel compiler is detected." )
  set( IntelCompiler_FOUND TRUE )

  get_dependency_version(IntelCompiler)

  if( "${IntelCompiler_VERSION}" LESS "1100" )
    message("Intel Compiler version ${IntelCompiler_VERSION} is not supported by CGAL-${CGAL_MAJOR_VERSION}.${CGAL_MINOR_VERSION}! (too old)")
  else()
    message( STATUS "Using Intel Compiler version 11 or later. Adding -fp-model strict" )
    uniquely_add_flags( CGAL_CXX_FLAGS "-fp-model strict" )
  endif()
endif()


if ( CMAKE_COMPILER_IS_GNUCXX )

  set( GCC_FOUND TRUE )
  
  get_dependency_version(GCC)
  
  if ( "${GCC_VERSION}" MATCHES "Not" OR "${GCC_VERSION}" MATCHES "Unknown" )
    set( GCC_FOUND FALSE )
  endif()
  
  if ( GCC_FOUND )
  
    if ( RUNNING_CGAL_AUTO_TEST )
      uniquely_add_flags( CGAL_CXX_FLAGS "-Wall" )
    endif()
    
    if ( "${GCC_VERSION}" MATCHES "^4." )
      message( STATUS "Using gcc version 4 or later. Adding -frounding-math" )
      uniquely_add_flags( CGAL_CXX_FLAGS "-frounding-math" )
    endif()
    
    if ( "${GCC_VERSION}" MATCHES "^4.2" )
      message( STATUS "Using gcc version 4.2. Adding -fno-strict-aliasing" )
      uniquely_add_flags( CGAL_CXX_FLAGS "-fno-strict-aliasing" )
    endif()
    
    if ( "${CMAKE_SYSTEM_PROCESSOR}" MATCHES "alpha" )
      message( STATUS "Using gcc on alpha. Adding -mieee -mfp-rounding-mode=d" )
      uniquely_add_flags( CGAL_CXX_FLAGS "-mieee -mfp-rounding-mode=d" )
    endif()
    
  endif()
  
endif()

include(CGAL_SetupFlags)

#--------------------------------------------------------------------------------------------------
#
#                                    -= DEPENDENCIES =-
#
#--------------------------------------------------------------------------------------------------

#
# The following variables are in the cache just so subdirectories can set them persistently.
# But they are not intended to persist from run to run as normal cache variables.
#
cache_set(CGAL_3RD_PARTY_DEFINITIONS    "" )
cache_set(CGAL_3RD_PARTY_INCLUDE_DIRS   "" )
cache_set(CGAL_3RD_PARTY_LIBRARIES      "" )
cache_set(CGAL_3RD_PARTY_LIBRARIES_DIRS "" ) 

cache_set(CGAL_Core_3RD_PARTY_DEFINITIONS    "" )
cache_set(CGAL_Core_3RD_PARTY_INCLUDE_DIRS   "" )
cache_set(CGAL_Core_3RD_PARTY_LIBRARIES      "" )
cache_set(CGAL_Core_3RD_PARTY_LIBRARIES_DIRS "" ) 

cache_set(CGAL_ImageIO_3RD_PARTY_DEFINITIONS    "" )
cache_set(CGAL_ImageIO_3RD_PARTY_INCLUDE_DIRS   "" )
cache_set(CGAL_ImageIO_3RD_PARTY_LIBRARIES      "" )
cache_set(CGAL_ImageIO_3RD_PARTY_LIBRARIES_DIRS "" )

cache_set(CGAL_Qt3_3RD_PARTY_DEFINITIONS    "" )
cache_set(CGAL_Qt3_3RD_PARTY_INCLUDE_DIRS   "" )
cache_set(CGAL_Qt3_3RD_PARTY_LIBRARIES      "" )
cache_set(CGAL_Qt3_3RD_PARTY_LIBRARIES_DIRS "" ) 

cache_set(CGAL_Qt4_3RD_PARTY_DEFINITIONS    "" )
cache_set(CGAL_Qt4_3RD_PARTY_INCLUDE_DIRS   "" )
cache_set(CGAL_Qt4_3RD_PARTY_LIBRARIES      "" )
cache_set(CGAL_Qt4_3RD_PARTY_LIBRARIES_DIRS "" ) 

cache_set(CGAL_PDB_3RD_PARTY_DEFINITIONS    "" )
cache_set(CGAL_PDB_3RD_PARTY_INCLUDE_DIRS   "" )
cache_set(CGAL_PDB_3RD_PARTY_LIBRARIES      "" )
cache_set(CGAL_PDB_3RD_PARTY_LIBRARIES_DIRS "" ) 

include(CGAL_SetupDependencies)

#--------------------------------------------------------------------------------------------------
#
#                                    -= Generation of compiler_config.h =-
#
#--------------------------------------------------------------------------------------------------

# The variables set are the #defines expected by compiler_config.h
# Note: CMake will not notice when files are added or removed
# but this is probably OK for the installation procedure.

file(GLOB all_config_tests "${CMAKE_SOURCE_DIR}/config/testfiles/*.cpp")

list( SORT all_config_tests )

file( WRITE ${CMAKE_BINARY_DIR}/include/CGAL/compiler_config.h "//\n// compiler_config.h is included by CGAL headers to load the needed compiler settings.\n//\n// DO NOT EDIT compiler_config.h. It is generated by CMake.\n//\n\n")

macro( add_config_flag flag )

    if ( ${flag} )
      file( APPEND ${CMAKE_BINARY_DIR}/include/CGAL/compiler_config.h "#define ${flag} 1\n\n")
    else()
      file( APPEND ${CMAKE_BINARY_DIR}/include/CGAL/compiler_config.h "//#define ${flag} 1\n\n")
    endif()

endmacro()

foreach(config_test_cpp ${all_config_tests})
    
    # Test's name is .cpp's base name
    get_filename_component(config_test_name ${config_test_cpp} NAME_WE)
    
    # Compile and run ${config_test_cpp}. Exit code is stored in ${config_test_name}.
    if( RUNNING_CGAL_AUTO_TEST OR "${config_test_name}" MATCHES "^${config_test_name}$" )

      CHECK_CXX_FILE_RUNS(${config_test_cpp} ${config_test_name} ${config_test_name})

      if ( ${config_test_name} )
        cache_set ( ${config_test_name} 0 )
      else()
        cache_set ( ${config_test_name} 1 )
      endif()  
      
    endif()
    
    add_config_flag( ${config_test_name} ${config_test_name} )
    
endforeach()

add_config_flag( CGAL_USE_GMP )
add_config_flag( CGAL_USE_MPFR )
add_config_flag( CGAL_USE_GMPXX )
add_config_flag( CGAL_USE_LEDA )

if ( RUNNING_CGAL_AUTO_TEST AND MSVC )
  file( APPEND "${CMAKE_BINARY_DIR}/include/CGAL/compiler_config.h" "#include <CGAL/Testsuite/vc_debug_hook.h>\n\n" )
endif()

#--------------------------------------------------------------------------------------------------
#
#                                    -= Installation Setup =-
#
#--------------------------------------------------------------------------------------------------

set ( CGAL_INSTALL_INC_DIR "include" CACHE STRING "The folder where CGAL header files will be installed, relative to CMAKE_INSTALL_PREFIX" )
set ( CGAL_INSTALL_LIB_DIR "lib"     CACHE STRING "The folder where CGAL libraries will be installed, relative to CMAKE_INSTALL_PREFIX" )

if ( NOT DEFINED CGAL_SOURCE_INSTALL )
  if ( WIN32 AND NOT UNIX AND NOT CYGWIN AND NOT CGAL_WIN32_CMAKE_ON_CYGWIN )
    set ( CGAL_SOURCE_INSTALL TRUE )
  else()
    set ( CGAL_SOURCE_INSTALL FALSE )
  endif()  
endif()

if ( CGAL_WIN32_CMAKE_ON_CYGWIN )
  exec_program(cygpath ARGS -w ${CMAKE_INSTALL_PREFIX} OUTPUT_VARIABLE CMAKE_INSTALL_PREFIX2 )
  file ( TO_CMAKE_PATH ${CMAKE_INSTALL_PREFIX2} CMAKE_INSTALL_PREFIX )
endif()

message( STATUS "CMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}" )

if ( CGAL_SOURCE_INSTALL )

  set ( CGAL_INSTALL_BIN_DIR "script"  CACHE STRING "The folder where CGAL user-side scripts will be installed, relative to CMAKE_INSTALL_PREFIX" )
  
  set ( CGAL_INSTALL_CMAKE_DIR "${CGAL_MODULES_REL_DIR}"  
        CACHE STRING "The folder where CGAL CMake modules will be installed, relative to CMAKE_INSTALL_PREFIX" 
      )
      
  set ( CGAL_INSTALL_ROOT_DIR "."  
        CACHE STRING "The folder where CGAL documentation and license files will be installed, relative to CMAKE_INSTALL_PREFIX" 
      )
      
  set ( CGAL_INSTALL_DOC_DIR "${CGAL_INSTALL_ROOT_DIR}" ) 
  
else()

  set ( CGAL_INSTALL_BIN_DIR    "bin"                     
        CACHE STRING "The folder where CGAL user-side scripts will be installed, relative to CMAKE_INSTALL_PREFIX" 
      )
      
  set ( CGAL_INSTALL_CMAKE_DIR "${CGAL_INSTALL_LIB_DIR}/CGAL"
        CACHE STRING "The folder where CGAL CMake modules will be installed, relative to CMAKE_INSTALL_PREFIX" 
      )
      
  set ( CGAL_INSTALL_DOC_DIR "share/doc/${CGAL_VERSION_DIR}"  
        CACHE STRING "The folder where CGAL documentation and license files will be installed, relative to CMAKE_INSTALL_PREFIX" 
      )

  set ( CGAL_INSTALL_MAN_DIR "share/man/man1"  
        CACHE STRING "The folder where manual pages for CGAL scripts will be installed, relative to CMAKE_INSTALL_PREFIX" 
      )

endif()

#--------------------------------------------------------------------------------------------------
#
#                                    -= Build =-
#
#--------------------------------------------------------------------------------------------------

set(CGAL_LIBRARIES_DIR ${CMAKE_BINARY_DIR}/lib)

set(LIBRARY_OUTPUT_PATH ${CGAL_LIBRARIES_DIR})

set(CGAL_INCLUDE_DIRS ${CMAKE_BINARY_DIR}/include ${CMAKE_SOURCE_DIR}/include)

include_directories (${CGAL_INCLUDE_DIRS})


add_subdirectory(src)

cache_get(CGAL_LIBRARY)
cache_get(CGAL_Core_LIBRARY)
cache_get(CGAL_ImageIO_LIBRARY)
cache_get(CGAL_PDB_LIBRARY)
cache_get(CGAL_Qt3_LIBRARY)
cache_get(CGAL_Qt4_LIBRARY)

cache_get(CGAL_3RD_PARTY_DEFINITIONS   )
cache_get(CGAL_3RD_PARTY_INCLUDE_DIRS  )
cache_get(CGAL_3RD_PARTY_LIBRARIES     )
cache_get(CGAL_3RD_PARTY_LIBRARIES_DIRS) 

cache_get(CGAL_Core_3RD_PARTY_DEFINITIONS   )
cache_get(CGAL_Core_3RD_PARTY_INCLUDE_DIRS  )
cache_get(CGAL_Core_3RD_PARTY_LIBRARIES     )
cache_get(CGAL_Core_3RD_PARTY_LIBRARIES_DIRS) 

cache_get(CGAL_ImageIO_3RD_PARTY_DEFINITIONS   )
cache_get(CGAL_ImageIO_3RD_PARTY_INCLUDE_DIRS  )
cache_get(CGAL_ImageIO_3RD_PARTY_LIBRARIES     )
cache_get(CGAL_ImageIO_3RD_PARTY_LIBRARIES_DIRS) 

cache_get(CGAL_Qt3_3RD_PARTY_DEFINITIONS   )
cache_get(CGAL_Qt3_3RD_PARTY_INCLUDE_DIRS  )
cache_get(CGAL_Qt3_3RD_PARTY_LIBRARIES     )
cache_get(CGAL_Qt3_3RD_PARTY_LIBRARIES_DIRS) 

cache_get(CGAL_Qt4_3RD_PARTY_DEFINITIONS   )
cache_get(CGAL_Qt4_3RD_PARTY_INCLUDE_DIRS  )
cache_get(CGAL_Qt4_3RD_PARTY_LIBRARIES     )
cache_get(CGAL_Qt4_3RD_PARTY_LIBRARIES_DIRS) 

cache_get(CGAL_PDB_3RD_PARTY_DEFINITIONS   )
cache_get(CGAL_PDB_3RD_PARTY_INCLUDE_DIRS  )
cache_get(CGAL_PDB_3RD_PARTY_LIBRARIES     )
cache_get(CGAL_PDB_3RD_PARTY_LIBRARIES_DIRS) 

hide_variable(EXECUTABLE_OUTPUT_PATH)
hide_variable(LIBRARY_OUTPUT_PATH)

#--------------------------------------------------------------------------------------------------
#
#                                    -= USER SIDE SETTINGS =-
#
#--------------------------------------------------------------------------------------------------

# Set CGAL_LIBRARY_NAME, CGAL_Qt3_LIBRARY_NAME and so on. Those variables
# are the name of CGAL libraries, without the path. Used in the generation
# of the installed CGALConfig.cmake
foreach(lib ${ALL_CGAL_LIBRARIES})
  get_filename_component(${lib}_LIBRARY_NAME "${${lib}_LIBRARY}" NAME CACHE)
  hide_variable(${lib}_LIBRARY_NAME)
endforeach()

create_CGALconfig_files()

#--------------------------------------------------------------------------------------------------
#
#                                    -= Installation Commands =-
#
#--------------------------------------------------------------------------------------------------

# WARNING: Use only relative paths; full paths break CPack!
# DESTINATION option is mandatory; skipping it breaks CPack!

install(FILES AUTHORS CHANGES LICENSE LICENSE.FREE_USE LICENSE.LGPL LICENSE.QPL DESTINATION ${CGAL_INSTALL_DOC_DIR} )


install(DIRECTORY include/CGAL                       DESTINATION ${CGAL_INSTALL_INC_DIR} )
install(DIRECTORY "${CMAKE_BINARY_DIR}/include/CGAL" DESTINATION ${CGAL_INSTALL_INC_DIR} )

file(GLOB scripts "scripts/*")
install(PROGRAMS ${scripts} DESTINATION ${CGAL_INSTALL_BIN_DIR})

install(DIRECTORY ${CGAL_MODULES_REL_DIR}/                DESTINATION ${CGAL_INSTALL_CMAKE_DIR} )
install(FILES     ${CGAL_MODULES_REL_DIR}/UseCGAL.cmake   DESTINATION ${CGAL_INSTALL_CMAKE_DIR} )


if ( GMP_IN_AUXILIARY )
  install(DIRECTORY auxiliary/gmp/include/ DESTINATION ${CGAL_INSTALL_INC_DIR} )
  install(DIRECTORY auxiliary/gmp/lib/     DESTINATION ${CGAL_INSTALL_LIB_DIR} )
endif()

if ( TAUCS_IN_AUXILIARY )
  install(DIRECTORY auxiliary/taucs/include/ DESTINATION ${CGAL_INSTALL_INC_DIR} )
  install(DIRECTORY auxiliary/tacus/lib/     DESTINATION ${CGAL_INSTALL_LIB_DIR} )
endif()

if ( ZLIB_IN_AUXILIARY )
  install(DIRECTORY auxiliary/zlib/include/ DESTINATION ${CGAL_INSTALL_INC_DIR} )
  install(DIRECTORY auxiliary/zlib/lib/     DESTINATION ${CGAL_INSTALL_LIB_DIR} )
endif()


if ( CGAL_SOURCE_INSTALL )
  install(FILES ${CMAKE_BINARY_DIR}/config/CGALConfig.cmake DESTINATION ${CGAL_INSTALL_ROOT_DIR} )
  install(FILES     CMakeLists.txt                          DESTINATION ${CGAL_INSTALL_ROOT_DIR} )
  install(DIRECTORY auxiliary                               DESTINATION ${CGAL_INSTALL_ROOT_DIR} )
  install(DIRECTORY src                                     DESTINATION ${CGAL_INSTALL_ROOT_DIR} )
else()
  install(FILES ${CMAKE_BINARY_DIR}/config/CGALConfig.cmake DESTINATION ${CGAL_INSTALL_CMAKE_DIR} )
  install(FILES auxiliary/cgal_create_cmake_script.1 DESTINATION ${CGAL_INSTALL_MAN_DIR} )
endif()


#--------------------------------------------------------------------------------------------------
#
#              -= Special installation commands to copy FindCGAL to the cmake directory =-
#
#--------------------------------------------------------------------------------------------------

# Installing FindCGAL into cmake itself is implemented as a custom target
# because it requires superuser privileges

if ( NOT "${CMAKE_ROOT}" STREQUAL "" )
  if ( EXISTS ${CMAKE_ROOT}/Modules )
    add_custom_target( install_FindCGAL
                      COMMAND ${CMAKE_COMMAND} -E copy ${CGAL_MODULES_DIR}/FindCGAL.cmake ${CMAKE_ROOT}/Modules
                    ) 
  endif()
endif()  

#--------------------------------------------------------------------------------------------------
#
#                                    -= APPLICATIONS =-
#
#--------------------------------------------------------------------------------------------------

macro( add_programs subdir target ON_OFF )

  cache_set( CGAL_EXECUTABLE_TARGETS "" )
  
  add_custom_target( ${target} )
  
  option( WITH_${target} "Select ${target}" ${ON_OFF} )
  if ( WITH_${target} )
    add_subdirectory( ${subdir} EXCLUDE_FROM_ALL )
  endif()
  
  cache_get( CGAL_EXECUTABLE_TARGETS )
  
  foreach( CGAL_EXECUTABLE_TARGET ${CGAL_EXECUTABLE_TARGETS} )
    add_dependencies( ${target} ${CGAL_EXECUTABLE_TARGET} )
  endforeach()
     
endmacro()

# This allows programs to locate CGALConfig.cmake
set(CGAL_DIR ${CMAKE_BINARY_DIR} )

if ( NOT RUNNING_CGAL_AUTO_TEST )

  add_programs(examples examples OFF )
  add_programs(demo     demos    OFF )
  
endif()

#--------------------------------------------------------------------------------------------------
#
#                                            -= CPack =-
#
#--------------------------------------------------------------------------------------------------

option( WITH_CPACK "Create package generation rules")
if( WITH_CPACK AND EXISTS "${CMAKE_ROOT}/Modules/CPack.cmake" )

    set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "CGAL - Computational Geometry Algorithms Library")
    set(CPACK_PACKAGE_VENDOR "CGAL Open Source Project")
    set(CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_SOURCE_DIR}/README")
    set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_SOURCE_DIR}/LICENSE")
    set(CPACK_PACKAGE_VERSION_MAJOR "${CGAL_MAJOR_VERSION}")
    set(CPACK_PACKAGE_VERSION_MINOR "${CGAL_MINOR_VERSION}")
    set(CPACK_PACKAGE_VERSION_PATCH "${CGAL_BUILD_VERSION}")
    set(CPACK_PACKAGE_INSTALL_DIRECTORY "CGAL ${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}")
    set(CPACK_SOURCE_PACKAGE_FILE_NAME "CGAL-${CGAL_VERSION}")
    set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_SOURCE_DIR}/LICENSE")
    
    if(NOT DEFINED CPACK_SYSTEM_NAME)
      set(CPACK_SYSTEM_NAME ${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR})
    endif()
    
    if(${CPACK_SYSTEM_NAME} MATCHES Windows)
      if(CMAKE_CL_64)
        set(CPACK_SYSTEM_NAME win64-${CMAKE_SYSTEM_PROCESSOR})
      else()
        set(CPACK_SYSTEM_NAME win32-${CMAKE_SYSTEM_PROCESSOR})
      endif()
    endif()
    
    if(NOT DEFINED CPACK_PACKAGE_FILE_NAME)
      set(CPACK_PACKAGE_FILE_NAME "${CPACK_SOURCE_PACKAGE_FILE_NAME}-${CPACK_SYSTEM_NAME}")
    endif()
    
    set(CPACK_PACKAGE_EXECUTABLES "CGAL" "CGAL")

    if(WIN32 AND NOT UNIX)
        set(CPACK_GENERATOR "NSIS")
        # There is a bug in NSI that does not handle full unix paths properly. Make
        # sure there is at least one set of four (4) backlasshes.
        #set(CPACK_PACKAGE_ICON "${CMAKE_SOURCE_DIR}\\\\cgal_install.gif")
        #set(CPACK_NSIS_INSTALLED_ICON_NAME "bin\\\\CGAL.exe")
        set(CPACK_NSIS_DISPLAY_NAME "${CPACK_PACKAGE_INSTALL_DIRECTORY} Computational Geometry Algorithms Library")
        set(CPACK_NSIS_HELP_LINK "http:\\\\\\\\www.cgal.org")
        set(CPACK_NSIS_URL_INFO_ABOUT "http:\\\\\\\\www.cgal.org")
        set(CPACK_NSIS_CONTACT "info@cgal.org")
        set(CPACK_NSIS_MODifY_PATH ON)
    else()
        set(CPACK_STRIP_FILES "")
        set(CPACK_SOURCE_STRIP_FILES "")
    endif()

    INCLUDE(CPack)
    
endif()

