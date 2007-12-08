# Top level CMakeLists.txt for CGAL

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

# The name of our project is "CGAL".  CMakeLists files in this project can
# refer to the root source directory of the project as ${CMAKE_SOURCE_DIR} or
# ${CGAL_SOURCE_DIR} and to the root binary directory of the project as
# ${CMAKE_BINARY_DIR} or ${CGAL_BINARY_DIR}.
project (CGAL)

if (NOT CMAKE_ROOT)
  message( FATAL_ERROR "CMAKE_ROOT enviroment variable not set")
endif()

# Check that the version of CMake is high enough.
# CPack was introduced in cmake 2.4
# FindQt3 is buggy in CMake 2.4.4.
CMAKE_MINIMUM_REQUIRED(VERSION 2.4.5 FATAL_ERROR)

# CGAL version number
set(CGAL_MAJOR_VERSION 3)
set(CGAL_MINOR_VERSION 4)
set(CGAL_BUILD_VERSION 1)
set(CGAL_VERSION
    "${CGAL_MAJOR_VERSION}.${CGAL_MINOR_VERSION}.${CGAL_BUILD_VERSION}")

# Just for fun
set(CMAKE_COLORMAKEFILE ON)


# # Debug: verbose output
# set(CMAKE_VERBOSE_MAKEFILE ON)
# set(CMAKE_START_TEMP_FILE "")
# set(CMAKE_END_TEMP_FILE "")

#-----------------------------------------------------------------------------
# TODO: build static and/or dynamic CGAL libraries?
# CGAL 3.2 builds only static libraries on Windows and static + dynamic on Unix.
if (WIN32)
    set(AUTO_LINK_ENABLED TRUE)
    set(BUILD_SHARED_LIBS OFF)
endif (WIN32)

#-----------------------------------------------------------------------------
# Configure install locations.
# WARNING: Use only relative paths (to ${CMAKE_INSTALL_PREFIX}); full paths break CPack!
set(CGAL_INSTALL_DIR		.)
set(CGAL_LIB_INSTALL_DIR	lib)
set(CGAL_INCLUDE_INSTALL_DIR	include)
set(CGAL_BIN_INSTALL_DIR	bin)
set(CGAL_AUXILIARY_INSTALL_DIR	auxiliary)
set(CGAL_MODULE_INSTALL_DIR	cmake/modules) # or ${CMAKE_ROOT}/Modules

#-----------------------------------------------------------------------------
# Where to look first for cmake modules, before ${CMAKE_ROOT}/Modules/ is checked
set(CMAKE_MODULE_PATH ${CGAL_SOURCE_DIR} ${CGAL_SOURCE_DIR}/cmake/modules)

# include CMake macros (without .cmake extension)
include(MacroOptionalFindPackage)
include(CheckCXXFileRuns)
include(LogCMakeVariables)
include(Qt3Macros)
include(AddFileDependencies)

#-----------------------------------------------------------------------------
# Attempt to use third party libraries (Boost, X11, GMP, GMPXX, MPFR, CORE,
# Zlib, LIDIA, LEDA, LEDAWIN, Qt, TAUCS and OpenGL)
# Create #define CGAL_USE_<library name> in compiler_config.h if found.
#
# Note: macro_optional_find_package() is the same as find_package(<name>) but additionally
# creates an OPTION(WITH_<name>) so the checking for the software can be disabled
# via ccmake or -DWITH_<name>=OFF

# Settings to compile an another project with CGAL
set(CGAL_3RD_PARTY_DEFINITIONS)
set(CGAL_3RD_PARTY_INCLUDE_DIRS)
set(CGAL_3RD_PARTY_LIBRARIES)

# Find Boost includes and optionaly Boost libraries
find_package(Boost2 REQUIRED)
if(Boost_FOUND)
    # Boost prints no message when found
    message(STATUS "Found Boost. Headers:  ${Boost_INCLUDE_DIRS}")
    message(STATUS "Found Boost. Libaries: ${Boost_LIBRARY_DIRS}")

    set(CGAL_3RD_PARTY_INCLUDE_DIRS ${CGAL_3RD_PARTY_INCLUDE_DIRS} ${Boost_INCLUDE_DIRS})
    set(CGAL_USE_BOOST 1)
endif(Boost_FOUND)

# Find BOOST_PROGRAM_OPTIONS
macro_optional_find_package(BOOST_PROGRAM_OPTIONS)
if(BOOST_PROGRAM_OPTIONS_FOUND)
    set(CGAL_3RD_PARTY_LIBRARIES ${CGAL_3RD_PARTY_LIBRARIES} ${BOOST_PROGRAM_OPTIONS_LIBRARIES})
    set(CGAL_USE_BOOST_PROGRAM_OPTIONS 1)
endif(BOOST_PROGRAM_OPTIONS_FOUND)

# Find X11 on Unix
macro_optional_find_package(X11)
if(X11_FOUND)
    # X11 prints no message when found
    message(STATUS "Found X11: ${X11_X11_LIB}")

    set(CGAL_3RD_PARTY_INCLUDE_DIRS ${CGAL_3RD_PARTY_INCLUDE_DIRS} ${X11_INCLUDE_DIR})
    set(CGAL_3RD_PARTY_LIBRARIES ${CGAL_3RD_PARTY_LIBRARIES} ${X11_X11_LIB})
    set(CGAL_USE_X11 1)
endif(X11_FOUND)

# Find GMP
macro_optional_find_package(GMP)
if(GMP_FOUND)
    set(CGAL_3RD_PARTY_INCLUDE_DIRS ${CGAL_3RD_PARTY_INCLUDE_DIRS} ${GMP_INCLUDE_DIR})
    set(CGAL_3RD_PARTY_LIBRARIES ${CGAL_3RD_PARTY_LIBRARIES} ${GMP_LIBRARIES})
    set(CGAL_USE_GMP 1)
endif(GMP_FOUND)

# Find GMPXX
macro_optional_find_package(GMPXX)
if(GMPXX_FOUND)
    set(CGAL_3RD_PARTY_INCLUDE_DIRS ${CGAL_3RD_PARTY_INCLUDE_DIRS} ${GMPXX_INCLUDE_DIR})
    set(CGAL_3RD_PARTY_LIBRARIES ${CGAL_3RD_PARTY_LIBRARIES} ${GMPXX_LIBRARIES})
    set(CGAL_USE_GMPXX 1)
endif(GMPXX_FOUND)

# Find MPFR
macro_optional_find_package(MPFR)
if(MPFR_FOUND)
    set(CGAL_3RD_PARTY_INCLUDE_DIRS ${CGAL_3RD_PARTY_INCLUDE_DIRS} ${MPFR_INCLUDE_DIR})
    set(CGAL_3RD_PARTY_LIBRARIES ${CGAL_3RD_PARTY_LIBRARIES} ${MPFR_LIBRARIES})
    set(CGAL_USE_MPFR 1)
endif(MPFR_FOUND)

# Find CORE in the system
macro_optional_find_package(CORE QUIET)
if(CORE_FOUND)
    message(STATUS "Found CORE: ${CORE_LIBRARIES}")

    set(CGAL_3RD_PARTY_INCLUDE_DIRS ${CGAL_3RD_PARTY_INCLUDE_DIRS} ${CORE_INCLUDE_DIR})
    set(CGAL_3RD_PARTY_LIBRARIES ${CGAL_3RD_PARTY_LIBRARIES} ${CORE_LIBRARIES})
    set(CGAL_USE_CORE 1)
else(CORE_FOUND)
    # If CORE not available in the system, use CORE shipped with CGAL
    macro_optional_find_package(CGAL_CORE)
    if(CGAL_CORE_FOUND)
	# Compile CORE shipped with CGAL
	add_subdirectory(src/CGALCore)

        set(CGAL_3RD_PARTY_INCLUDE_DIRS ${CGAL_3RD_PARTY_INCLUDE_DIRS} ${CGAL_CORE_INCLUDE_DIR})
        set(CGAL_3RD_PARTY_LIBRARIES ${CGAL_3RD_PARTY_LIBRARIES} ${CGAL_CORE_LIBRARIES})
        set(CGAL_USE_CGAL_CORE 1)
        set(CGAL_USE_CORE 1)

        # by compatibility with install_cgal
        set(CGAL_3RD_PARTY_DEFINITIONS ${CGAL_3RD_PARTY_DEFINITIONS} "-DCGAL_USE_CGAL_CORE")
    endif(CGAL_CORE_FOUND)
endif(CORE_FOUND)

# Find Zlib
macro_optional_find_package(ZLIB)
if(ZLIB_FOUND)
    set(CGAL_3RD_PARTY_INCLUDE_DIRS ${CGAL_3RD_PARTY_INCLUDE_DIRS} ${ZLIB_INCLUDE_DIR})
    set(CGAL_3RD_PARTY_LIBRARIES ${CGAL_3RD_PARTY_LIBRARIES} ${ZLIB_LIBRARIES})
    set(CGAL_USE_ZLIB 1)
endif(ZLIB_FOUND)

# TODO: Write FindLIDIA.cmake

# TODO: Write FindLEDA.cmake

# TODO: Write FindLEDAWIN.cmake

# Find Qt3
macro_optional_find_package(Qt3)
if(QT_FOUND)
    # Qt3 prints no message when found
    message(STATUS "Found Qt3: ${QT_QT_LIBRARY}")

    set(CGAL_3RD_PARTY_INCLUDE_DIRS ${CGAL_3RD_PARTY_INCLUDE_DIRS} ${QT_INCLUDE_DIR})
    set(CGAL_3RD_PARTY_LIBRARIES ${CGAL_3RD_PARTY_LIBRARIES} ${QT_QT_LIBRARY})
    set(CGAL_3RD_PARTY_DEFINITIONS ${CGAL_3RD_PARTY_DEFINITIONS} ${QT_DEFINITIONS})
    set(CGAL_USE_QT 1)
endif(QT_FOUND)

# Find TAUCS
macro_optional_find_package(TAUCS)
if(TAUCS_FOUND)
    set(CGAL_3RD_PARTY_INCLUDE_DIRS ${CGAL_3RD_PARTY_INCLUDE_DIRS} ${TAUCS_INCLUDE_DIR})
    set(CGAL_3RD_PARTY_LIBRARIES ${CGAL_3RD_PARTY_LIBRARIES} ${TAUCS_LIBRARIES})
    set(CGAL_USE_TAUCS 1)
endif()

# Find OpenGL
macro_optional_find_package(OpenGL)
if(OPENGL_FOUND)
    # OpenGL prints no message when found
    message(STATUS "Found OpenGL: ${OPENGL_LIBRARIES}")

    set(CGAL_3RD_PARTY_INCLUDE_DIRS ${CGAL_3RD_PARTY_INCLUDE_DIRS} ${OPENGL_INCLUDE_DIR})
    set(CGAL_3RD_PARTY_LIBRARIES ${CGAL_3RD_PARTY_LIBRARIES} ${OPENGL_LIBRARIES})
    set(CGAL_USE_OpenGL 1) # unused by install_cgal
endif(OPENGL_FOUND)

#-----------------------------------------------------------------------------
# Installation.

# Install files from top level folder
# WARNING: Use only relative paths; full paths break CPack!
#          DESTINATION option is mandatory; skipping it breaks CPack!
INSTALL(FILES CHANGES
	      INSTALL INSTALL.MacOSX INSTALL.win32
	      LICENSE LICENSE.FREE_USE LICENSE.LGPL LICENSE.QPL
	      README
	      VERSION
	DESTINATION ${CGAL_INSTALL_DIR})

# Install headers
add_subdirectory(include/CGAL)

# Install scripts
add_subdirectory(scripts)

# Install CMake modules
add_subdirectory(cmake/modules)

# Install precompiled third-party libraries (Windows only).
add_subdirectory(auxiliary)

# #-----------------------------------------------------------------------------
# # Save the compiler settings so another project can import them.
# CMAKE_EXPORT_BUILD_SETTINGS(${CGAL_BINARY_DIR}/CGALBuildSettings.cmake)
# INSTALL_FILES(${CGAL_MODULE_INSTALL_DIR} .cmake CGALBuildSettings)

#-----------------------------------------------------------------------------
# Platform configuration tests.
# The variables set are the #defines expected by compiler_config.h

# For each .cpp configuration test
# Note: CMake will not notice when files are added or removed
# but this is probably OK for the installation procedure.
FILE(GLOB all_config_tests "${CGAL_SOURCE_DIR}/config/testfiles/*.cpp")
FOREACH(config_test_cpp ${all_config_tests})
    # Test's name is .cpp's base name
    GET_FILENAME_COMPONENT(config_test_name ${config_test_cpp} NAME_WE)

    # Compile and run ${config_test_cpp}. Exit code is stored in ${config_test_name}.
    CHECK_CXX_FILE_RUNS(${config_test_cpp} ${config_test_name} ${config_test_name})

    # compiler_config.h should #define ${config_test_name} is the test failed
    # => flip ${config_test_name}
    if (${config_test_name})
      set(${config_test_name} 0)
    else()
      set(${config_test_name} 1)
    endif()
ENDFOREACH(config_test_cpp ${all_config_tests})

#-----------------------------------------------------------------------------
# TODO: Rewrite here install_cgal's create_standard_header_replacements()

#-----------------------------------------------------------------------------
# Save third-party libraries' and configuration tests' #defines for use
# when building CGAL libraries.
CONFIGURE_FILE(${CGAL_SOURCE_DIR}/compiler_config.h.in
               ${CGAL_BINARY_DIR}/include/CGAL/compiler_config.h @ONLY IMMEDIATE)
INSTALL(FILES ${CGAL_BINARY_DIR}/include/CGAL/compiler_config.h DESTINATION ${CGAL_INCLUDE_INSTALL_DIR}/CGAL)

message(FATAL_ERROR "stop")

#-----------------------------------------------------------------------------
# CGAL requires special compiler flags on some platforms.
# 1) set settings to compile CGAL libraries.

# Execute install_cgal's set_compiler_flags().
include(set_compiler_flags)

# Include directories
set(CGAL_INCLUDE_DIRS ${CGAL_BINARY_DIR}/include
                      ${CGAL_SOURCE_DIR}/include
                      ${CGAL_3RD_PARTY_INCLUDE_DIRS})

# Compiler settings
set(CGAL_DEFINITIONS ${CGAL_3RD_PARTY_DEFINITIONS} ${ADDITIONAL_CXXFLAGS})

# Set Microsoft Visual C++ compilation warning level
if (CMAKE_CXX_COMPILER MATCHES "^(CL|cl)")
    # CGAL 3.2 level is 2
    SET(CMAKE_CXX_WARNING_LEVEL 2
	CACHE STRING "C++ compiler warning level" FORCE)
    MARK_AS_ADVANCED(CMAKE_CXX_WARNING_LEVEL)
endif (CMAKE_CXX_COMPILER MATCHES "^(CL|cl)")

# Settings to compile CGAL libraries
if (BUILD_SHARED_LIBS)
    SET(CGAL_LIB_CXX_FLAGS ${CGAL_DEFINITIONS} ${CGAL_SHARED_LIB_CXXFLAGS})
else (BUILD_SHARED_LIBS)
    SET(CGAL_LIB_CXX_FLAGS ${CGAL_DEFINITIONS} ${CGAL_LIB_CXXFLAGS})
endif (BUILD_SHARED_LIBS)

# Settings to link CGAL libraries
SET(CGAL_LIB_EXE_LINKER_FLAGS "${ADDITIONAL_LDFLAGS} ${CGAL_LIB_LDFLAGS}")
SET(CGAL_LIB_SHARED_LINKER_FLAGS "${ADDITIONAL_LDFLAGS} ${CGAL_SHARED_LIB_LDFLAGS}")
if (MATH_LIB)
    set(CGAL_3RD_PARTY_LIBRARIES ${CGAL_3RD_PARTY_LIBRARIES} "-l${MATH_LIB}")
endif (MATH_LIB)

#-----------------------------------------------------------------------------
# Build and install CGAL and CGALQt libraries

# Folder where the compiler will put CGAL libraries
set(LIBRARY_OUTPUT_PATH ${CGAL_BINARY_DIR}/lib)

# Append the library version information to the library target properties.
set(CGAL_LIBRARY_PROPERTIES ${CGAL_LIBRARY_PROPERTIES}
    VERSION "${CGAL_VERSION}"
    SOVERSION "${CGAL_MAJOR_VERSION}.${CGAL_MINOR_VERSION}"
    )

# Build and install CGAL library
add_subdirectory(src/CGAL)

if(QT_FOUND)
    # Build and install CGALQt library
    add_subdirectory(src/CGALQt)
endif(QT_FOUND)

#-----------------------------------------------------------------------------
# CGAL requires special compiler flags on some platforms.
# 2) set settings to use (i.e. compile with) CGAL libraries.

# Make sure the compiler can find include files
include_directories (${CGAL_INCLUDE_DIRS})

# Add compiler flags CGAL needs to work on this platform.
ADD_DEFINITIONS(${CGAL_DEFINITIONS})

# Settings to link with CGAL libraries
#
# common linker settings
set(CGAL_LIBRARIES
    ${ADDITIONAL_LDFLAGS}
)
# add CGAL library
GET_TARGET_PROPERTY(CGAL_LOCATION CGAL LOCATION)
set(CGAL_LIBRARIES ${CGAL_LIBRARIES}
    ${CGAL_LOCATION}
)
# add CGALQt library
if(QT_FOUND)
    GET_TARGET_PROPERTY(CGALQt_LOCATION CGALQt LOCATION)
    set(CGAL_LIBRARIES ${CGAL_LIBRARIES}
	${CGALQt_LOCATION}
    )
endif(QT_FOUND)
# add third-party libraries
set(CGAL_LIBRARIES ${CGAL_LIBRARIES}
    ${CGAL_3RD_PARTY_LIBRARIES}
)

#-----------------------------------------------------------------------------
# TODO: Configure Dart testing support

set(BUILD_TESTING ON CACHE BOOL "Turn testing on/off")
# INCLUDE(Dart)
# MARK_AS_ADVANCED(DART_ROOT TCL_TCLSH)
if(BUILD_TESTING)
    ENABLE_TESTING()
#     CONFIGURE_FILE(${CGAL_SOURCE_DIR}/CMake/CTestCustom.ctest.in
#         ${CGAL_BINARY_DIR}/CMake/CTestCustom.ctest @ONLY)
#     FILE(WRITE ${CGAL_BINARY_DIR}/CTestCustom.cmake
#         "INCLUDE(\"${CGAL_BINARY_DIR}/CMake/CTestCustom.ctest\")\n")
endif(BUILD_TESTING)

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
            endif(_file MATCHES "CMakeLists\\.txt")
        ENDFOREACH(_file)
    ENDFOREACH(_folder demo examples test)
endif(BUILD_TESTING)

# #-----------------------------------------------------------------------------
# # Help other projects use CGAL.
#
# # write the compiler settings and the install dirs
# # into a cmake file which will be installed
# CONFIGURE_FILE(${CGAL_SOURCE_DIR}/UseCGAL.cmake.in
#                ${CGAL_BINARY_DIR}/UseCGAL.cmake @ONLY IMMEDIATE)
# INSTALL_FILES(${CGAL_MODULE_INSTALL_DIR} .cmake UseCGAL)

#-----------------------------------------------------------------------------
# Build a CPack installer
# CPack was introduced in cmake 2.4
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

#-----------------------------------------------------------------------------
#if(WIN32)
#  MESSAGE("Configure done. You may need to click on configure again if the OK button is disabled.")
#endif(WIN32)

#-----------------------------------------------------------------------------
# DEBUG
#LOG_CMAKE_VARIABLES("end")
