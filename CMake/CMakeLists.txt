# Top level CMakeLists.txt for CGAL
# The name of our project is "CGAL".  CMakeLists files in this project can
# refer to the root source directory of the project as ${CMAKE_SOURCE_DIR} or
# ${CGAL_SOURCE_DIR} and to the root binary directory of the project as
# ${CMAKE_BINARY_DIR} or ${CGAL_BINARY_DIR}.
project (CGAL CXX)

#--------------------------------------------------------------------------------------------------
#
#                                    -= BASIC SETUP =-
#
#--------------------------------------------------------------------------------------------------

# Where to look first for cmake modules, before ${CMAKE_ROOT}/Modules/ is checked
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CGAL_SOURCE_DIR}/cmake/modules)

include(CGALcommon)
include(GeneratorSpecificSettings)
include(CheckCXXFileRuns)
#include(Qt3Macros)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING
      "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel."
      FORCE
     )
else()
  set(CMAKE_CONFIGURATION_TYPES ${CMAKE_BUILD_TYPE} CACHE STRING
      "The types of build for multiple-configuration generators"
      FORCE
     )
endif()




# CGAL version number
# TODO: This should be parsed from the VERSION file
set(CGAL_MAJOR_VERSION 3)
set(CGAL_MINOR_VERSION 4)
set(CGAL_BUILD_VERSION 1)
set(CGAL_VERSION "${CGAL_MAJOR_VERSION}.${CGAL_MINOR_VERSION}.${CGAL_BUILD_VERSION}")


# WARNING: Use only relative paths; full paths break CPack!
set(CGAL_INSTALL_DIR		        .)
set(CGAL_LIB_INSTALL_DIR	      lib)
set(CGAL_INCLUDE_INSTALL_DIR	  include)
set(CGAL_AUXILIARY_INSTALL_DIR	auxiliary)
set(CGAL_CMAKE_INSTALL_DIR     	cmake/modules)



#--------------------------------------------------------------------------------------------------
#
#                                    -= DEPENDENCIES =-
#
#--------------------------------------------------------------------------------------------------

# Settings to compile an another project with CGAL
set(CGAL_3RD_PARTY_DEFINITIONS)
set(CGAL_3RD_PARTY_INCLUDE_DIRS)
set(CGAL_3RD_PARTY_LIBRARIES)
set(CGAL_3RD_PARTY_LIBRARIES_DIR)

# This module finds dependent libraries: Boost, GMP, CORE, etc..
# Updates CGAL_3RD_PARTY_* accordingly
# Sets CGAL_USE_* which is stored in compiler_config.h
# Each dependency can be selected or unselected via the option WITH_*
include(FindCGALDependencies)



#--------------------------------------------------------------------------------------------------
#
#                                    -= Generation of compiler_config.h =-
#
#--------------------------------------------------------------------------------------------------

# Loads the contents of the configuration file to check that it contains all needed flags
file( READ "${CGAL_SOURCE_DIR}/compiler_config.h.in" COMPILER_CONFIG )

# The variables set are the #defines expected by compiler_config.h
# For each .cpp configuration test
# Note: CMake will not notice when files are added or removed
# but this is probably OK for the installation procedure.

file(GLOB all_config_tests "${CGAL_SOURCE_DIR}/config/testfiles/*.cpp")

foreach(config_test_cpp ${all_config_tests})
    
    # Test's name is .cpp's base name
    get_filename_component(config_test_name ${config_test_cpp} NAME_WE)
    
    if ( NOT "${COMPILER_CONFIG}" MATCHES "${config_test_name}" )
      message( SEND_ERROR "File 'compiler_config.h.in' is not up to date. Flag ${config_test_name} is missing!" )
    endif()

    # Compile and run ${config_test_cpp}. Exit code is stored in ${config_test_name}.
    CHECK_CXX_FILE_RUNS(${config_test_cpp} ${config_test_name} ${config_test_name})

    # compiler_config.h should #define ${config_test_name} is the test failed
    # => flip ${config_test_name}
    if (${config_test_name})
      set(${config_test_name} 0)
    else()
      set(${config_test_name} 1)
    endif()
endforeach()

configure_file(${CGAL_SOURCE_DIR}/compiler_config.h.in ${CGAL_SOURCE_DIR}/include/CGAL/compiler_config.h @ONLY IMMEDIATE)
               

#--------------------------------------------------------------------------------------------------
#
#                                    -= Build =-
#
#--------------------------------------------------------------------------------------------------



set(CGAL_INCLUDE_DIRS ${CGAL_SOURCE_DIR}/include)

set(CGAL_DEFINITIONS ${COMMON_DEFINITIONS} )

SET(CGAL_LIB_SHARED_LINKER_FLAGS "${ADDITIONAL_LDFLAGS} ${CGAL_SHARED_LIB_LDFLAGS}")

add_subdirectory(src/CGAL)

if( CGAL_CORE_FOUND )
  add_subdirectory(src/CGALCore)
endif()

if(QT_FOUND)
  add_subdirectory(src/CGALQt)
endif()

add_subdirectory(examples/Straight_skeleton_2)


#-----------------------------------------------------------------------------
# TODO: Configure Dart testing support

set(BUILD_TESTING OFF CACHE BOOL "Turn testing on/off")

# INCLUDE(Dart)
# MARK_AS_ADVANCED(DART_ROOT TCL_TCLSH)
if(BUILD_TESTING)
    ENABLE_TESTING()
#     configure_file(${CGAL_SOURCE_DIR}/CMake/CTestCustom.ctest.in ${CGAL_BINARY_DIR}/CMake/CTestCustom.ctest @ONLY)
#     FILE(WRITE ${CGAL_BINARY_DIR}/CTestCustom.cmake "INCLUDE(\"${CGAL_BINARY_DIR}/CMake/CTestCustom.ctest\")\n")
endif()

#-----------------------------------------------------------------------------
# Recurse into demo/example/test subdirectories
# TODO: move compilation of demos/examples/tests from "make all" to "make test"

if(BUILD_TESTING)
  # Find all folders under demo/example/test which contain a CMakeLists.txt.
  # TODO: rewrite this section because CMake will not notice when packages
  # are added or removed.
  FOREACH(_folder demo examples test)
    FILE(GLOB_RECURSE all_files "${_folder}/*")
    FOREACH(_file ${all_files})
      if(_file MATCHES "CMakeLists\\.txt")
        GET_FILENAME_COMPONENT(_path ${_file} PATH)
        #MESSAGE(STATUS "Configure '${_path}'")
        add_subdirectory(${_path})
      endif()
    ENDFOREACH()
  ENDFOREACH()
endif()


#--------------------------------------------------------------------------------------------------
#
#                                    -= INSTALLATION =-
#
#--------------------------------------------------------------------------------------------------

# DESTINATION option is mandatory; skipping it breaks CPack!

install(FILES CHANGES "INSTALL" INSTALL.MacOSX INSTALL.win32
	            LICENSE LICENSE.FREE_USE LICENSE.LGPL LICENSE.QPL
	            README "VERSION"  
	      DESTINATION ${CGAL_INSTALL_DIR}
       )

install(DIRECTORY "include/CGAL" DESTINATION ${CGAL_INCLUDE_INSTALL_DIR} )
install(DIRECTORY "auxiliary"    DESTINATION ${CGAL_AUXILIARY_INSTALL_DIR} )
install(DIRECTORY "cmake"        DESTINATION ${CGAL_CMAKE_INSTALL_DIR} )




#--------------------------------------------------------------------------------------------------
#
#                                    -= USER SIDE SETTINGS =-
#
#--------------------------------------------------------------------------------------------------

# FinCGAL is paltform independet
configure_file(${CGAL_SOURCE_DIR}/FindCGAL.cmake.in ${CGAL_SOURCE_DIR}/cmake/modules/FindCGAL.cmake @ONLY IMMEDIATE)

# UseCGAL OTOH contains platform depedent definitions
configure_file(${CGAL_SOURCE_DIR}/UseCGAL.cmake.in  ${CGAL_BINARY_DIR}/cmake/modules/UseCGAL.cmake  @ONLY IMMEDIATE)

install( FILES ${CGAL_BINARY_DIR}/cmake/modules/FindCGAL.cmake ${CGAL_BINARY_DIR}/cmake/modules/FindCGAL.cmake
         DESTINATION ${CMAKE_ROOT}/Modules
       )  


#--------------------------------------------------------------------------------------------------
#
#                                            -= CPack =-
#
#--------------------------------------------------------------------------------------------------

IF(EXISTS "${CMAKE_ROOT}/Modules/CPack.cmake")
    SET(CPACK_PACKAGE_DESCRIPTION_SUMMARY "CGAL - Computational Geometry Algorithms Library")
    SET(CPACK_PACKAGE_VENDOR "CGAL Open Source Project")
    SET(CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_SOURCE_DIR}/README")
    SET(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_SOURCE_DIR}/LICENSE")
    SET(CPACK_PACKAGE_VERSION_MAJOR "${CGAL_MAJOR_VERSION}")
    SET(CPACK_PACKAGE_VERSION_MINOR "${CGAL_MINOR_VERSION}")
    SET(CPACK_PACKAGE_VERSION_PATCH "${CGAL_BUILD_VERSION}")
    SET(CPACK_PACKAGE_INSTALL_DIRECTORY "CGAL ${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}")
    SET(CPACK_SOURCE_PACKAGE_FILE_NAME "CGAL-${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}")
    SET(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_SOURCE_DIR}/LICENSE")
    IF(NOT DEFINED CPACK_SYSTEM_NAME)
    SET(CPACK_SYSTEM_NAME ${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR})
    ENDIF(NOT DEFINED CPACK_SYSTEM_NAME)
    IF(${CPACK_SYSTEM_NAME} MATCHES Windows)
    IF(CMAKE_CL_64)
      SET(CPACK_SYSTEM_NAME win64-${CMAKE_SYSTEM_PROCESSOR})
    ELSE(CMAKE_CL_64)
      SET(CPACK_SYSTEM_NAME win32-${CMAKE_SYSTEM_PROCESSOR})
    ENDIF(CMAKE_CL_64)
    ENDIF(${CPACK_SYSTEM_NAME} MATCHES Windows)
    IF(NOT DEFINED CPACK_PACKAGE_FILE_NAME)
      SET(CPACK_PACKAGE_FILE_NAME "${CPACK_SOURCE_PACKAGE_FILE_NAME}-${CPACK_SYSTEM_NAME}")
    ENDIF(NOT DEFINED CPACK_PACKAGE_FILE_NAME)
    #SET(CPACK_PACKAGE_EXECUTABLES "CGAL" "CGAL")

    IF(WIN32 AND NOT UNIX)
        SET(CPACK_GENERATOR "NSIS")
        # There is a bug in NSI that does not handle full unix paths properly. Make
        # sure there is at least one set of four (4) backlasshes.
        SET(CPACK_PACKAGE_ICON "${CMAKE_SOURCE_DIR}\\\\cgal_install.gif")
        #SET(CPACK_NSIS_INSTALLED_ICON_NAME "bin\\\\CGAL.exe")
        SET(CPACK_NSIS_DISPLAY_NAME "${CPACK_PACKAGE_INSTALL_DIRECTORY} Computational Geometry Algorithms Library")
        SET(CPACK_NSIS_HELP_LINK "http:\\\\\\\\www.cgal.org")
        SET(CPACK_NSIS_URL_INFO_ABOUT "http:\\\\\\\\www.cgal.com")
        SET(CPACK_NSIS_CONTACT "info@cgal.org")
        SET(CPACK_NSIS_MODIFY_PATH ON)
    ELSE(WIN32 AND NOT UNIX)
        SET(CPACK_STRIP_FILES "")
        SET(CPACK_SOURCE_STRIP_FILES "")
    ENDIF(WIN32 AND NOT UNIX)

    INCLUDE(CPack)
ENDIF(EXISTS "${CMAKE_ROOT}/Modules/CPack.cmake")

