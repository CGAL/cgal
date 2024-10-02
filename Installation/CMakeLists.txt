# Top level CMakeLists.txt for CGAL
# The name of our project is "CGAL".  CMakeLists files in this project can
# refer to the root source directory of the project as ${CGAL_SOURCE_DIR} or
# ${CGAL_SOURCE_DIR} and to the root binary directory of the project as
# ${CGAL_BINARY_DIR} or ${CGAL_BINARY_DIR}.
if(NOT PROJECT_NAME)
  cmake_minimum_required(VERSION 3.12...3.29)
  project(CGAL CXX C)
endif()

# Fix a bug: `try_compile` should use `CMAKE_EXE_LINKER_FLAGS`.
# That variable can contain important flags like
#
#     -static-libstdc++ -static-libgcc
#
# Use GNUInstallDirs to get canonical paths
include(GNUInstallDirs)
#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------
#
#                                    -= PACKAGE SETUP =-
#
#--------------------------------------------------------------------------------------------------

if(CGAL_BRANCH_BUILD)

  message(STATUS "Build CGAL from ${CGAL_SCM_NAME}-branch: ${CGAL_SCM_BRANCH_NAME}")

  # list packages
  file(
    GLOB CGAL_CONFIGURED_PACKAGES
    RELATIVE ${CGAL_SOURCE_DIR}
    "${CGAL_SOURCE_DIR}/*")
  list(SORT CGAL_CONFIGURED_PACKAGES)
  list(REMOVE_ITEM CGAL_CONFIGURED_PACKAGES copyright CMakeLists.txt .svn .git)

  # detect and remove not existing package-directories
  foreach(package ${CGAL_CONFIGURED_PACKAGES})
    if(NOT IS_DIRECTORY "${CGAL_SOURCE_DIR}/${package}")
      list(APPEND CGAL_WRONG_PACKAGES ${package})
    elseif(NOT EXISTS
           "${CGAL_SOURCE_DIR}/${package}/package_info/${package}/maintainer")
      list(APPEND CGAL_WRONG_PACKAGES ${package})
    endif()
  endforeach()

  find_program(
    CGAL_CREATE_CMAKE_SCRIPT cgal_create_cmake_script HINT
    ${CGAL_SOURCE_DIR}/Scripts/scripts
    DOC "Script cgal_create_cmake_script, that creates CMakeLists.txt files"
    NO_DEFAULT_PATH NO_CMAKE_ENVIRONMENT_PATH)

  if(DEFINED CGAL_WRONG_PACKAGES)
    message(STATUS "Removed not-a-package: ${CGAL_WRONG_PACKAGES}")
    list(REMOVE_ITEM CGAL_CONFIGURED_PACKAGES ${CGAL_WRONG_PACKAGES})
  endif()

  set(CGAL_CONFIGURED_PACKAGES_NAMES ${CGAL_CONFIGURED_PACKAGES})
  set(CGAL_CONFIGURED_PACKAGES)
  foreach(package ${CGAL_CONFIGURED_PACKAGES_NAMES})
    list(APPEND CGAL_CONFIGURED_PACKAGES "${CGAL_SOURCE_DIR}/${package}")
  endforeach()

  set(CGAL_INSTALLATION_PACKAGE_DIR
      "${CGAL_SOURCE_DIR}/Installation"
      CACHE INTERNAL "Directory containing the Installation package")
  set(CGAL_MAINTENANCE_PACKAGE_DIR
      "${CGAL_SOURCE_DIR}/Maintenance"
      CACHE INTERNAL "Directory containing the Maintenance package")
  set(CGAL_CORE_PACKAGE_DIR
      "${CGAL_SOURCE_DIR}/Core"
      CACHE INTERNAL "Directory containing the Core package")
  set(CGAL_GRAPHICSVIEW_PACKAGE_DIR
      "${CGAL_SOURCE_DIR}/GraphicsView"
      CACHE INTERNAL "Directory containing the GraphicsView package")
  file(GLOB INCLUDE_CGAL_QT_DIRECTORIES "${CGAL_SOURCE_DIR}/*/include/CGAL/Qt")
  message(STATUS "Include CGAL/Qt directories: ${INCLUDE_CGAL_QT_DIRECTORIES}")

  message(STATUS "Installation package directory: ${CGAL_INSTALLATION_PACKAGE_DIR}")
  message(STATUS "Maintenance package directory: ${CGAL_MAINTENANCE_PACKAGE_DIR}")
  message(STATUS "Core package directory: ${CGAL_CORE_PACKAGE_DIR}")

else(CGAL_BRANCH_BUILD)

  # get release name
  file(TO_CMAKE_PATH ${CGAL_SOURCE_DIR} CMAKE_SOURCE_CDIR)
  string(REGEX REPLACE "^/" "" CMAKE_SOURCE_DDIR ${CMAKE_SOURCE_CDIR})
  string(REPLACE "/" ";" CMAKE_SOURCE_PDIR ${CMAKE_SOURCE_DDIR})
  list(GET CMAKE_SOURCE_PDIR -1 CGAL_RELEASE_NAME)
  message(STATUS "Build CGAL from release in directory ${CGAL_RELEASE_NAME}")

  set(CGAL_CONFIGURED_PACKAGES_NAMES ${CGAL_RELEASE_NAME})

  # list packages
  set(CGAL_CONFIGURED_PACKAGES "${CGAL_SOURCE_DIR}")

  set(CGAL_INSTALLATION_PACKAGE_DIR
      "${CGAL_SOURCE_DIR}"
      CACHE INTERNAL "Directory containing the Installation package")
  set(CGAL_MAINTENANCE_PACKAGE_DIR
      "${CGAL_SOURCE_DIR}"
      CACHE INTERNAL "Directory containing the Maintenance package")
  set(CGAL_CORE_PACKAGE_DIR
      "${CGAL_SOURCE_DIR}"
      CACHE INTERNAL "Directory containing the Core package")
  set(CGAL_GRAPHICSVIEW_PACKAGE_DIR
      "${CGAL_SOURCE_DIR}"
      CACHE INTERNAL "Directory containing the GraphicsView package")
  set(INCLUDE_CGAL_QT_DIRECTORIES "${CGAL_SOURCE_DIR}/include/CGAL/Qt")

# Enable testing with CGAL_ENABLE_TESTING. Before CGAL-6.0, users would enable
# the tests by specifying BUILD_TESTING. For compatibility, If BUILD_TESTING is
# set, that is the default value for CGAL_ENABLE_TESTING. Otherwise, the default
# value is OFF.
option(CGAL_ENABLE_TESTING "Build the testing tree." ${BUILD_TESTING})
  if(CGAL_ENABLE_TESTING)
    enable_testing()
  endif()
endif(CGAL_BRANCH_BUILD)

#allow to force disabling boost multiprecision support
option(CGAL_DO_NOT_USE_BOOST_MP "Disable the support of boost multiprecision library" FALSE)

#message(STATUS "Packages found: ${CGAL_CONFIGURED_PACKAGES}")

list(SORT CGAL_CONFIGURED_PACKAGES_NAMES)
if(IS_DIRECTORY "${CGAL_SOURCE_DIR}/Documentation")
  # Rescope CGAL_CONFIGURED_PACKAGES_NAMES
  # We need this variable in the Doxygen configuration process. Reset it
  # and change the scope to the parent scope.
  # Philipp Moeller, 2013-05-302013-05-302013-05-30
  set(CGAL_CONFIGURED_PACKAGES_NAMES
      ${CGAL_CONFIGURED_PACKAGES_NAMES}
      PARENT_SCOPE)
endif()

message(STATUS "Packagenames: ${CGAL_CONFIGURED_PACKAGES_NAMES}")
message("== Setting paths (DONE) ==\n")

message("== Generate version files ==")

if(CGAL_BRANCH_BUILD)

  #
  # Create version files
  #

  file(REMOVE ${CGAL_BINARY_DIR}/VERSION)
  if(CGAL_BUGFIX_VERSION AND CGAL_BUGFIX_VERSION GREATER 0)
    set(CGAL_CREATED_VERSION_NUM
        "${CGAL_MAJOR_VERSION}.${CGAL_MINOR_VERSION}.${CGAL_BUGFIX_VERSION}")
  else()
    set(CGAL_CREATED_VERSION_NUM "${CGAL_MAJOR_VERSION}.${CGAL_MINOR_VERSION}")
  endif()

  # Accessed in the documentation.
  set(CGAL_CREATED_VERSION_NUM
      ${CGAL_CREATED_VERSION_NUM}
      PARENT_SCOPE)

  # TODO EBEB set number of internal release (replace "900")
  set(CGAL_BUILD_VERSION 900)
  set(CGAL_CREATED_VERSION
      "${CGAL_CREATED_VERSION_NUM}-I-${CGAL_BUILD_VERSION}")
  set(CGAL_VERSION "${CGAL_CREATED_VERSION_NUM}-I-${CGAL_BUILD_VERSION}")

  file(WRITE "${CGAL_BINARY_DIR}/VERSION" "${CGAL_CREATED_VERSION_NUM}")

  message(STATUS "CGAL_VERSION is ${CGAL_CREATED_VERSION}")
  set(CGAL_FULL_VERSION
      ${CGAL_CREATED_VERSION}
      CACHE INTERNAL "CGAL version")

  #define CGAL_VERSION_NR 1030500135
  string(LENGTH "${CGAL_MAJOR_VERSION}" CGAL_MAJOR_LEN)
  if("${CGAL_MAJOR_LEN}" EQUAL 1)
    set(CGAL_TMP_MAJOR "0${CGAL_MAJOR_VERSION}")
  endif()
  string(LENGTH "${CGAL_MINOR_VERSION}" CGAL_MINOR_LEN)
  if("${CGAL_MINOR_LEN}" EQUAL 1)
    set(CGAL_TMP_MINOR "0${CGAL_MINOR_VERSION}")
  else()
    set(CGAL_TMP_MINOR "${CGAL_MINOR_VERSION}")
  endif()

  set(CGAL_TMP_BUILD "${CGAL_BUILD_VERSION}")
  string(LENGTH "${CGAL_TMP_BUILD}" CGAL_BUILD_LEN)

  while(CGAL_BUILD_LEN LESS 4)
    set(CGAL_TMP_BUILD "0${CGAL_TMP_BUILD}")
    string(LENGTH "${CGAL_TMP_BUILD}" CGAL_BUILD_LEN)
  endwhile()

  set(CGAL_CREATED_VERSION_NR
      "1${CGAL_TMP_MAJOR}${CGAL_TMP_MINOR}${CGAL_BUGFIX_VERSION}${CGAL_TMP_BUILD}"
  )
  message(STATUS "CGAL_VERSION_NR is ${CGAL_CREATED_VERSION_NR}")

  message(STATUS "CGAL_GIT_HASH is ${CGAL_GIT_HASH}")
  message(
    STATUS "CGAL_CREATED_SVN_REVISION is ${CGAL_CREATED_SVN_REVISION} (dummy)")

  file(REMOVE ${CGAL_BINARY_DIR}/include/CGAL/version.h)
  string(TIMESTAMP CGAL_CREATED_RELEASE_DATE "%Y%m%d")
  configure_file(
    ${CGAL_INSTALLATION_PACKAGE_DIR}/cmake/modules/config/version.h.in
    ${CGAL_BINARY_DIR}/include/CGAL/version.h @ONLY)

  #
  # Duplicate files
  #

  if(CGAL_REPORT_DUPLICATE_FILES)

    message("== Searching for duplicate files ==")

    foreach(package ${CGAL_CONFIGURED_PACKAGES})

      file(GLOB_RECURSE CGAL_PACKAGE_HEADER_FILES FOLLOW_SYMLINKS
           ${package}/*.h)
      list(SORT CGAL_PACKAGE_HEADER_FILES)
      foreach(file ${CGAL_PACKAGE_HEADER_FILES})
        string(REPLACE "${package}/" "" pure_file ${file})
        if(NOT ${pure_file} MATCHES ".*CMakeFiles.*")
          list(APPEND CGAL_HEADER_FILES ${pure_file})
        endif()
      endforeach()

      file(GLOB_RECURSE CGAL_PACKAGE_CPP_FILES FOLLOW_SYMLINKS ${package}/*.cpp)
      list(SORT CGAL_PACKAGE_CPP_FILES)
      foreach(file ${CGAL_PACKAGE_CPP_FILES})
        string(REPLACE "${package}/" "" pure_file ${file})
        if(NOT ${pure_file} MATCHES ".*CMakeFiles.*")
          list(APPEND CGAL_CPP_FILES ${pure_file})
        endif()
      endforeach()

    endforeach()

    list(SORT CGAL_HEADER_FILES)
    set(CGAL_UNIQUE_HEADER_FILES ${CGAL_HEADER_FILES})
    set(CGAL_DUPLICATE_HEADER_FILES ${CGAL_HEADER_FILES})
    list(REMOVE_DUPLICATES CGAL_UNIQUE_HEADER_FILES)
    #message(STATUS "Headers: ${CGAL_HEADER_FILES}")
    foreach(file ${CGAL_UNIQUE_HEADER_FILES})
      list(FIND CGAL_DUPLICATE_HEADER_FILES "${file}" CGAL_FILE_IDX)
      list(REMOVE_AT CGAL_DUPLICATE_HEADER_FILES ${CGAL_FILE_IDX})
    endforeach()
    message(
      STATUS
        "Duplicate header files (will be repeated at end): ${CGAL_DUPLICATE_HEADER_FILES}"
    )

    list(SORT CGAL_CPP_FILES)
    set(CGAL_UNIQUE_CPP_FILES ${CGAL_CPP_FILES})
    set(CGAL_DUPLICATE_CPP_FILES ${CGAL_CPP_FILES})
    list(REMOVE_DUPLICATES CGAL_UNIQUE_CPP_FILES)
    #message(STATUS "CPPs: ${CGAL_CPP_FILES}")
    foreach(file ${CGAL_UNIQUE_CPP_FILES})
      list(FIND CGAL_DUPLICATE_CPP_FILES "${file}" CGAL_FILE_IDX)
      list(REMOVE_AT CGAL_DUPLICATE_CPP_FILES ${CGAL_FILE_IDX})
    endforeach()
    message(
      STATUS
        "Duplicate cpp files (will be repeated at end): ${CGAL_DUPLICATE_CPP_FILES}"
    )

  endif(CGAL_REPORT_DUPLICATE_FILES)

else()

  if(NOT EXISTS "${CGAL_SOURCE_DIR}/VERSION"
     OR NOT EXISTS "${CGAL_SOURCE_DIR}/include/CGAL/version.h")
    message(
      FATAL_ERROR
        "File VERSION or include/CGAL/version.h are missing, required by release-build."
    )
  endif()

  #
  # Read and parse CGAL version number from VERSION file
  #

  file(READ "${CGAL_SOURCE_DIR}/VERSION" CGAL_VERSION)
  # Comment from Laurent Rineau, 2012/06/28:
  #   file(READ ..) should be replace by file(STRINGS ..), and then we take
  #   the first line. That would avoid parsing errors when the VERSION file
  #   contains a CR character!

  set(CGAL_FULL_VERSION ${CGAL_VERSION})
  string(REPLACE "-" "." CGAL_VERSION_TOKENS1 ${CGAL_VERSION})
  string(REPLACE "." ";" CGAL_VERSION_TOKENS ${CGAL_VERSION_TOKENS1})
  list(LENGTH CGAL_VERSION_TOKENS CGAL_VERSION_TOKENS_LEN)
  list(GET CGAL_VERSION_TOKENS 0 CGAL_MAJOR_VERSION)

  set(CGAL_MINOR_VERSION)
  set(CGAL_BUGFIX_VERSION)
  set(CGAL_BUILD_VERSION 1000)

  if(CGAL_VERSION_TOKENS_LEN GREATER 1)
    list(GET CGAL_VERSION_TOKENS 1 CGAL_MINOR_VERSION)
  endif()

  if(CGAL_VERSION_TOKENS_LEN GREATER 2)
    list(GET CGAL_VERSION_TOKENS 2 CGAL_BUGFIX_VERSION)

    # The following condition will be true the token #2 is not a number (if
    # it is for example "I" or "Ic"):
    if(NOT CGAL_BUGFIX_VERSION GREATER 0 AND NOT CGAL_BUGFIX_VERSION EQUAL 0)
      set(CGAL_BUGFIX_VERSION 0)
      if(CGAL_VERSION_TOKENS_LEN GREATER 3)
        list(GET CGAL_VERSION_TOKENS 3 CGAL_BUILD_VERSION)
      endif()
    else()
      if(CGAL_VERSION_TOKENS_LEN GREATER 4)
        list(GET CGAL_VERSION_TOKENS 4 CGAL_BUILD_VERSION)
      endif()
    endif()
    # If CGAL_BUILD_VERSION is not a strictly positive number, error
    if(NOT CGAL_BUILD_VERSION GREATER 0)
      message(
        WARNING
          "Error parsing VERSION file: CGAL_BUILD_VERSION is not a number.\n"
          "The content of the VERSION file is:\n${CGAL_VERSION}\n"
          "CGAL_BUILD_VERSION was set to:\n${CGAL_BUILD_VERSION}")
      set(CGAL_BUILD_VERSION 1000)
    endif()
  endif()

  if(NOT CGAL_BUGFIX_VERSION)
    set(CGAL_BUGFIX_VERSION 0)
  endif()

  if(CGAL_BUGFIX_VERSION AND CGAL_BUGFIX_VERSION GREATER 0)
    set(CGAL_CREATED_VERSION_NUM
        "${CGAL_MAJOR_VERSION}.${CGAL_MINOR_VERSION}.${CGAL_BUGFIX_VERSION}")
  else()
    set(CGAL_CREATED_VERSION_NUM "${CGAL_MAJOR_VERSION}.${CGAL_MINOR_VERSION}")
  endif()
endif()

#--------------------------------------------------------------------------------------------------
#
#                                    -= BASIC SETUP =-
#
#--------------------------------------------------------------------------------------------------

message(STATUS "CGAL_MAJOR_VERSION=${CGAL_MAJOR_VERSION}")
message(STATUS "CGAL_MINOR_VERSION=${CGAL_MINOR_VERSION}")
message(STATUS "CGAL_BUGFIX_VERSION=${CGAL_BUGFIX_VERSION}")
if(CGAL_BUILD_VERSION LESS 1000)
  message(STATUS "CGAL_BUILD_VERSION=${CGAL_BUILD_VERSION}")
endif()

set(CGAL_VERSION_DIR CGAL-${CGAL_VERSION})
set(CGAL_MODULES_REL_DIR cmake/modules)
set(CGAL_MODULES_DIR ${CGAL_INSTALLATION_PACKAGE_DIR}/${CGAL_MODULES_REL_DIR})

include(${CGAL_MODULES_DIR}/CGAL_Macros.cmake)
include(${CGAL_MODULES_DIR}/CGAL_enable_end_of_configuration_hook.cmake)
cgal_setup_module_path()

if(RUNNING_CGAL_AUTO_TEST OR CGAL_TEST_SUITE)
  message(STATUS "Operating system:")
  execute_process(
    COMMAND uname -a
    TIMEOUT 5
    OUTPUT_VARIABLE uname_a
    ERROR_VARIABLE uname_a)
  message(STATUS "${uname_a}")
  cgal_display_compiler_version()
endif()

include(CGAL_Common)
include(CGAL_GeneratorSpecificSettings)
include(CGAL_CheckCXXFileRuns)
include(CGAL_add_test)

mark_as_advanced(CMAKE_BACKWARDS_COMPATIBILITY)

message(
  STATUS
    "USING CMake version: ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}.${CMAKE_PATCH_VERSION}"
)
message(STATUS "System: ${CMAKE_SYSTEM_NAME}")

#--------------------------------------------------------------------------------------------------
#
#                                    -= FLAGS =-
#
#--------------------------------------------------------------------------------------------------

if(MSVC)

  uniquely_add_flags(CGAL_CXX_FLAGS
                     "-D_SCL_SECURE_NO_DEPRECATE;-D_SCL_SECURE_NO_WARNINGS")
  uniquely_add_flags(CGAL_CXX_FLAGS "/fp:strict")
  uniquely_add_flags(CGAL_CXX_FLAGS "/fp:except-")
  uniquely_add_flags(CGAL_CXX_FLAGS "/bigobj") # Use /bigobj by default

  if(RUNNING_CGAL_AUTO_TEST OR CGAL_TEST_SUITE)
    set(CMAKE_CXX_WARNING_LEVEL
        2
        CACHE STRING "MSVC C++ compiler warning level" FORCE)
    mark_as_advanced(CMAKE_CXX_WARNING_LEVEL)
  endif()

endif()

if("\"${CMAKE_CXX_COMPILER_ID}\"" MATCHES "SunPro")
  message(STATUS "Using SunPro compiler, using STLPort 4.")

  uniquely_add_flags(CGAL_CXX_FLAGS
                     "-features=extensions;-library=stlport4;-D_GNU_SOURCE")
  uniquely_add_flags(CGAL_SHARED_LINKER_FLAGS "-library=stlport4")
  uniquely_add_flags(CGAL_EXE_LINKER_FLAGS "-library=stlport4")
endif()

if("${CMAKE_CXX_COMPILER}" MATCHES "icl" OR "${CMAKE_CXX_COMPILER}" MATCHES
                                            "icpc")
  message(STATUS "Intel compiler is detected.")
  set(IntelCompiler_FOUND TRUE)

  get_dependency_version(IntelCompiler)

  if("${IntelCompiler_VERSION}" LESS "1100")
    message(
      "Intel Compiler version ${IntelCompiler_VERSION} is not supported by CGAL-${CGAL_MAJOR_VERSION}.${CGAL_MINOR_VERSION}! (too old, must be 11.0 or after)"
    )
  else()
    message(
      STATUS "Using Intel Compiler version 11 or later. Adding -fp-model strict"
    )
    if(WIN32)
      uniquely_add_flags(CGAL_CXX_FLAGS "/fp:strict")
    else()
      uniquely_add_flags(CGAL_CXX_FLAGS "-fp-model strict")
    endif()
  endif()
endif()

if(CMAKE_COMPILER_IS_GNUCXX)

  set(GCC_FOUND TRUE)

  get_dependency_version(GCC)

  if("${GCC_VERSION}" MATCHES "Not" OR "${GCC_VERSION}" MATCHES "Unknown")
    set(GCC_FOUND FALSE)
  endif()

  if(GCC_FOUND)

    if(RUNNING_CGAL_AUTO_TEST OR CGAL_TEST_SUITE)
      uniquely_add_flags(CGAL_CXX_FLAGS "-Wall")
      # Remove -g from the relevant CMAKE_CXX_FLAGS. This will also
      # propagate to the rest of the tests, since we overwrite those
      # flags with the ones used to build CGAL.
      string(REGEX REPLACE "-g( |$)" "" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
      # We only allow the release types DEBUG and RELEASE, but handle
      # all possible values just to be sure.
      foreach(release_type DEBUG RELEASE MINSIZEREL RELWITHDEBINFO)
        string(REGEX REPLACE "-g( |$)" "" CMAKE_CXX_FLAGS_${release_type}
                             "${CMAKE_CXX_FLAGS_${release_type}}")
      endforeach()

    endif()

    if("${GCC_VERSION}" MATCHES "^[4-9].")
      message(STATUS "Using gcc version 4 or later. Adding -frounding-math")
      uniquely_add_flags(CGAL_CXX_FLAGS "-frounding-math")
    endif()

    if("${CMAKE_SYSTEM_PROCESSOR}" MATCHES "alpha")
      message(STATUS "Using gcc on alpha. Adding -mieee -mfp-rounding-mode=d")
      uniquely_add_flags(CGAL_CXX_FLAGS "-mieee -mfp-rounding-mode=d")
    endif()

  endif()

endif()

if(CGAL_TEST_SUITE)
  include(display-third-party-libs-versions)
endif()

message("== Generate version files (DONE) ==\n")
#--------------------------------------------------------------------------------------------------
#
#                                    -= External libraries =-
#
#--------------------------------------------------------------------------------------------------

if(CGAL_DEV_MODE OR RUNNING_CGAL_AUTO_TEST OR CGAL_TEST_SUITE)
  message("== Set up flags ==")

  # Ugly hack to be compatible with current CGAL testsuite process (as of
  # Nov. 2017). -- Laurent Rineau
  include(CGAL_SetupFlags)

  message("== Set up flags (DONE) ==\n")
else()
  include(CGAL_display_flags)
endif()

message("== Detect external libraries ==")

# set some to have special prefix
macro(set_special_prefix library prefix)
  set(SPECIAL_PREFIXES
      "${SPECIAL_PREFIXES}set(CGAL_EXT_LIB_${library}_PREFIX \"${prefix}\")\n")
  set(CGAL_EXT_LIB_${library}_PREFIX ${prefix})
endmacro()

set_special_prefix(Qt6 QT)
set_special_prefix(Eigen3 EIGEN3)
set_special_prefix(Coin3D COIN3D)

# some libraries are essential (stl and Boost.Thread are treated in another way)
if($ENV{CGAL_DISABLE_GMP})
  set(CGAL_DISABLE_GMP
      ON
      CACHE INTERNAL "")
endif()
if(CGAL_DISABLE_GMP)
  message("Disable the GMP support")
  list(LENGTH CGAL_ESSENTIAL_3RD_PARTY_LIBRARIES CGAL_ESSENTIAL_LENGTH)
  if(NOT CGAL_ESSENTIAL_LENGTH EQUAL 0)
    list(REMOVE_ITEM CGAL_ESSENTIAL_3RD_PARTY_LIBRARIES GMP MPFR)
  endif()
  unset(CGAL_CONFIGURED_LIBRARIES)
  unset(CGAL_CONFIGURED_LIBRARIES CACHE)
  unset(GMP_FOUND)
  unset(GMP_FOUND CACHE)
  unset(GMPXX_FOUND)
  unset(GMPXX_FOUND CACHE)
  unset(MPFR_FOUND)
  unset(MPFR_FOUND CACHE)
  unset(WITH_CGAL_Core)
  unset(WITH_CGAL_Core CACHE)
  unset(CGAL_USE_GMP)
  unset(CGAL_USE_GMP CACHE)
  # Nasty trick to make sure <gmp.h> is not used when CGAL_DISABLE_GMP is TRUE
  file(WRITE "${CGAL_BINARY_DIR}/include/gmp.h"
       "#error GMP is disabled by the CMake option CGAL_DISABLE_GMP")
else()
  list(APPEND CGAL_ESSENTIAL_3RD_PARTY_LIBRARIES GMP MPFR)

  # When CMake is run several times, to avoid duplicates
  list(REMOVE_DUPLICATES CGAL_ESSENTIAL_3RD_PARTY_LIBRARIES)

  file(REMOVE "${CGAL_BINARY_DIR}/include/gmp.h")
endif()
hide_variable(CGAL_ESSENTIAL_3RD_PARTY_LIBRARIES)

#--------------------------------------------------------------------------------------------------
#
#                                    -= DEPENDENCIES =-
#
#--------------------------------------------------------------------------------------------------

#
# The following variables are in the cache just so subdirectories can set them persistently.
# But they are not intended to persist from run to run as normal cache variables.
# Similar variables are created when a library is detected.
cache_set(CGAL_3RD_PARTY_DEFINITIONS "")
cache_set(CGAL_3RD_PARTY_INCLUDE_DIRS "")
cache_set(CGAL_3RD_PARTY_LIBRARIES "")
cache_set(CGAL_3RD_PARTY_LIBRARIES_DIRS "")

# additional info: some header files in CGAL add additional code if
# certain optional libs are installed, and some examples/tests rely on
# this; e.g. in MPFI/RS in Algebraic_kernel_d. For these cases CGAL
# and the example/test must be configured with MPFI (just one is not sufficient)

message("== Detect external libraries (DONE) ==\n")

#--------------------------------------------------------------------------------------------------
#
#                                    -= Generation of compiler_config.h =-
#
#--------------------------------------------------------------------------------------------------

#--------------------------------------------------------------------------------------------------
#
#                                    -= Installation Setup =-
#
#--------------------------------------------------------------------------------------------------

set(CGAL_INSTALL_INC_DIR
    "${CMAKE_INSTALL_INCLUDEDIR}"
    CACHE
      STRING
      "The folder where CGAL header files will be installed, relative to CMAKE_INSTALL_PREFIX"
)
set(CGAL_INSTALL_LIB_DIR
    "${CMAKE_INSTALL_LIBDIR}"
    CACHE
      STRING
      "The folder where CGAL libraries will be installed, relative to CMAKE_INSTALL_PREFIX"
)

if(CGAL_WIN32_CMAKE_ON_CYGWIN)
  execute_process(COMMAND cygpath -w ${CMAKE_INSTALL_PREFIX}
                  OUTPUT_STRIP_TRAILING_WHITESPACE
                  OUTPUT_VARIABLE CMAKE_INSTALL_PREFIX2)
  file(TO_CMAKE_PATH ${CMAKE_INSTALL_PREFIX2} CMAKE_INSTALL_PREFIX)
endif()

set(CGAL_INSTALL_BIN_DIR
    "${CMAKE_INSTALL_BINDIR}"
    CACHE
      STRING
      "The folder where CGAL user-side scripts will be installed, relative to CMAKE_INSTALL_PREFIX"
)

set(CGAL_INSTALL_CMAKE_DIR
    "${CGAL_INSTALL_LIB_DIR}/cmake/CGAL"
    CACHE
      STRING
      "The folder where CGAL CMake modules will be installed, relative to CMAKE_INSTALL_PREFIX"
)

set(CGAL_INSTALL_DOC_DIR
    "${CMAKE_INSTALL_DOCDIR}"
    CACHE
      STRING
      "The folder where CGAL documentation and license files will be installed, relative to CMAKE_INSTALL_PREFIX"
)

set(CGAL_INSTALL_MAN_DIR
    "${CMAKE_INSTALL_MANDIR}/man1"
    CACHE
      STRING
      "The folder where manual pages for CGAL scripts will be installed, relative to CMAKE_INSTALL_PREFIX"
)

#--------------------------------------------------------------------------------------------------
#
#                                    -= Build =-
#
#--------------------------------------------------------------------------------------------------

message("== Generating build files ==")

set(CGAL_INCLUDE_DIRS ${CGAL_BINARY_DIR}/include)

foreach(package ${CGAL_CONFIGURED_PACKAGES})
  if(EXISTS ${package}/include)
    set(CGAL_INCLUDE_DIRS ${CGAL_INCLUDE_DIRS} ${package}/include)
  endif()
endforeach()

include_directories(${CGAL_INCLUDE_DIRS})

cache_get(CGAL_3RD_PARTY_DEFINITIONS   )
cache_get(CGAL_3RD_PARTY_INCLUDE_DIRS  )
cache_get(CGAL_3RD_PARTY_LIBRARIES     )
cache_get(CGAL_3RD_PARTY_LIBRARIES_DIRS)

install(DIRECTORY ${INCLUDE_CGAL_QT_DIRECTORIES} DESTINATION "${CGAL_INSTALL_INC_DIR}/CGAL/" COMPONENT CGAL_Qt6)
if(CGAL_BRANCH_BUILD)
  install(DIRECTORY "${CGAL_GRAPHICSVIEW_PACKAGE_DIR}/demo/resources/" DESTINATION "${CGAL_INSTALL_CMAKE_DIR}/demo/resources" COMPONENT CGAL_Qt6)
  install(DIRECTORY "${CGAL_GRAPHICSVIEW_PACKAGE_DIR}/demo/icons/" DESTINATION "${CGAL_INSTALL_CMAKE_DIR}/demo/icons" COMPONENT CGAL_Qt6)
endif()

#
# Variables used when WITH_{demos|examples|tests} are TRUE
#
get_property(CGAL_FOUND GLOBAL PROPERTY CGAL_FOUND)
get_property(CGAL_Core_FOUND GLOBAL PROPERTY CGAL_Core_FOUND)
get_property(CGAL_ImageIO_FOUND GLOBAL PROPERTY CGAL_ImageIO_FOUND)
get_property(CGAL_Qt6_FOUND GLOBAL PROPERTY CGAL_Qt6_FOUND)

#
# Repeat some problems
#
message("== Generating build files (DONE) ==\n")

if(CGAL_BRANCH_BUILD AND CGAL_REPORT_DUPLICATE_FILES)

  message(STATUS "Problems: ")

  if(CGAL_DUPLICATE_HEADER_FILES)
    message(WARNING "WARNING: Duplicate header files")
    foreach(hfile ${CGAL_DUPLICATE_HEADER_FILES})
      message(STATUS "  File '${hfile}' multiply found in ")
      foreach(package ${CGAL_CONFIGURED_PACKAGES})
        file(GLOB_RECURSE CGAL_PACKAGE_HEADER_FILES FOLLOW_SYMLINKS
             ${package}/*.h)
        list(SORT CGAL_PACKAGE_HEADER_FILES)
        foreach(file ${CGAL_PACKAGE_HEADER_FILES})
          string(REPLACE "${package}/" "" pure_file ${file})
          if("${pure_file}" STREQUAL "${hfile}")
            message(STATUS "    ${package}")
          endif()
        endforeach()
      endforeach()
    endforeach()
  endif()

  if(CGAL_DUPLICATE_CPP_FILES)
    message(WARNING "WARNING: Duplicate cpp files")
    foreach(cfile ${CGAL_DUPLICATE_CPP_FILES})
      message(STATUS "  File '${cfile}' multiply found in ")
      foreach(package ${CGAL_CONFIGURED_PACKAGES})
        file(GLOB_RECURSE CGAL_PACKAGE_CPP_FILES FOLLOW_SYMLINKS
             ${package}/*.cpp)
        list(SORT CGAL_PACKAGE_CPP_FILES)
        foreach(file ${CGAL_PACKAGE_CPP_FILES})
          string(REPLACE "${package}/" "" pure_file ${file})
          if("${pure_file}" STREQUAL "${cfile}")
            message(STATUS "    ${package}")
          endif()
        endforeach()
      endforeach()
    endforeach()
  endif()

endif()

configure_file("${CGAL_MODULES_DIR}/CGALConfig_binary_header_only.cmake.in"
  "${CGAL_BINARY_DIR}/CGALConfig.cmake"  @ONLY)
configure_file("${CGAL_MODULES_DIR}/CGALConfigVersion_binary_header_only.cmake.in"
  "${CGAL_BINARY_DIR}/CGALConfigVersion.cmake"  @ONLY)

#--------------------------------------------------------------------------------------------------
#
#                                    -= Installation Commands =-
#
#--------------------------------------------------------------------------------------------------

# DESTINATION option is mandatory; skipping it breaks CPack!

if(CGAL_INSTALL_DOC_DIR)
  install(FILES AUTHORS CHANGES.md LICENSE LICENSE.GPL
                LICENSE.LGPL DESTINATION ${CGAL_INSTALL_DOC_DIR})
endif()

#install all includes collected in trunk et cetera

foreach(dir ${CGAL_CONFIGURED_PACKAGES})
  if(EXISTS ${dir}/include/CGAL)
    install(
      DIRECTORY ${dir}/include/CGAL
      DESTINATION ${CGAL_INSTALL_INC_DIR}
      PATTERN ".svn" EXCLUDE
      PATTERN "Qt" EXCLUDE)
  endif()
endforeach()
if(EXISTS ${CGAL_BINARY_DIR}/include/CGAL)
  install(
    DIRECTORY ${CGAL_BINARY_DIR}/include/CGAL
    DESTINATION ${CGAL_INSTALL_INC_DIR}
    PATTERN ".svn" EXCLUDE
    PATTERN "Qt" EXCLUDE)
endif()

file(GLOB scripts "scripts/*")
list(SORT scripts)
list(REMOVE_ITEM scripts
     ${CMAKE_CURRENT_SOURCE_DIR}/scripts/cgal_create_assertions.sh
     ${CMAKE_CURRENT_SOURCE_DIR}/scripts/cgal_create_cmake_script_with_options)
install(PROGRAMS ${scripts} DESTINATION ${CGAL_INSTALL_BIN_DIR})

install(DIRECTORY ${CGAL_MODULES_REL_DIR}/
        DESTINATION ${CGAL_INSTALL_CMAKE_DIR})

if(IS_DIRECTORY auxiliary/gmp/include AND IS_DIRECTORY auxiliary/gmp/lib)
  install(DIRECTORY auxiliary/gmp/include/ DESTINATION ${CGAL_INSTALL_INC_DIR})
  install(DIRECTORY auxiliary/gmp/lib/ DESTINATION ${CGAL_INSTALL_LIB_DIR})
endif()

if(ZLIB_IN_AUXILIARY)
  install(DIRECTORY auxiliary/zlib/include/ DESTINATION ${CGAL_INSTALL_INC_DIR})
  install(DIRECTORY auxiliary/zlib/lib/ DESTINATION ${CGAL_INSTALL_LIB_DIR})
endif()

install(FILES
  ${CMAKE_CURRENT_LIST_DIR}/lib/cmake/CGAL/CGALConfig.cmake
  ${CMAKE_CURRENT_LIST_DIR}/lib/cmake/CGAL/CGALConfigBuildVersion.cmake
  ${CMAKE_CURRENT_LIST_DIR}/lib/cmake/CGAL/CGALConfigVersion.cmake
  DESTINATION ${CGAL_INSTALL_CMAKE_DIR})
if(NOT CGAL_INSTALL_CMAKE_DIR STREQUAL "${CGAL_INSTALL_LIB_DIR}/cmake/CGAL")
  configure_file(${CMAKE_CURRENT_LIST_DIR}/lib/cmake/CGAL/CGALConfig-installation-dirs.cmake.in
    ${CGAL_BINARY_DIR}/config/CGALConfig-installation-dirs.cmake)
  install(FILES
    ${CGAL_BINARY_DIR}/config/CGALConfig-installation-dirs.cmake
    DESTINATION ${CGAL_INSTALL_CMAKE_DIR})
endif()

if(CGAL_INSTALL_MAN_DIR)
  install(FILES auxiliary/cgal_create_cmake_script.1
          DESTINATION ${CGAL_INSTALL_MAN_DIR})
endif()

#--------------------------------------------------------------------------------------------------
#
#              -= Special installation commands to copy FindCGAL to the cmake directory =-
#
#--------------------------------------------------------------------------------------------------

# Installing FindCGAL into cmake itself is implemented as a custom target
# because it requires superuser privileges

if(NOT "${CMAKE_ROOT}" STREQUAL "")
  if(EXISTS ${CMAKE_ROOT}/Modules)
    add_custom_target(
      install_FindCGAL
      COMMAND ${CMAKE_COMMAND} -E copy ${CGAL_MODULES_DIR}/FindCGAL.cmake
              ${CMAKE_ROOT}/Modules)
  endif()
endif()

#--------------------------------------------------------------------------------------------------
#
#                                    -= APPLICATIONS =-
#
#--------------------------------------------------------------------------------------------------

if(NOT TARGET ALL_CGAL_TARGETS)
  add_custom_target(ALL_CGAL_TARGETS)
endif()

macro(add_programs subdir target ON_OFF)

  cache_set(CGAL_EXECUTABLE_TARGETS "")

  add_custom_target(${target})
  add_dependencies(ALL_CGAL_TARGETS ${target})

  option(WITH_${target} "Select ${target}" ${ON_OFF})
  if(WITH_${target})
    add_subdirectory(${subdir} EXCLUDE_FROM_ALL)
  endif()

  cache_get(CGAL_EXECUTABLE_TARGETS)

  foreach(CGAL_EXECUTABLE_TARGET ${CGAL_EXECUTABLE_TARGETS})
    add_dependencies(${target} ${CGAL_EXECUTABLE_TARGET})
  endforeach()

endmacro()

# This allows programs to locate CGALConfig.cmake
set(CGAL_DIR ${CGAL_BINARY_DIR})

if(NOT RUNNING_CGAL_AUTO_TEST)

  add_programs(examples examples OFF)
  add_programs(demo demos OFF)
  if(IS_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/test")
    add_programs(test tests OFF)
  endif()
endif()

message("== Setting paths ==")

#--------------------------------------------------------------------------------------------------
#
#                                            -= CPack =-
#
#--------------------------------------------------------------------------------------------------

#option( WITH_CPACK "Create package generation rules")
if(WITH_CPACK AND EXISTS "${CMAKE_ROOT}/Modules/CPack.cmake")

  set(CPACK_PACKAGE_DESCRIPTION_SUMMARY
      "CGAL - Computational Geometry Algorithms Library")
  set(CPACK_PACKAGE_VENDOR "CGAL Open Source Project")
  set(CPACK_PACKAGE_DESCRIPTION_FILE "${CGAL_INSTALLATION_PACKAGE_DIR}/README")
  set(CPACK_RESOURCE_FILE_LICENSE "${CGAL_INSTALLATION_PACKAGE_DIR}/LICENSE")
  set(CPACK_PACKAGE_VERSION_MAJOR "${CGAL_MAJOR_VERSION}")
  set(CPACK_PACKAGE_VERSION_MINOR "${CGAL_MINOR_VERSION}")
  set(CPACK_PACKAGE_VERSION_PATCH "${CGAL_BUGFIX_VERSION}")
  set(CPACK_PACKAGE_INSTALL_DIRECTORY
      "CGAL ${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}")
  set(CPACK_SOURCE_PACKAGE_FILE_NAME "CGAL-${CGAL_VERSION}")
  set(CPACK_RESOURCE_FILE_LICENSE "${CGAL_INSTALLATION_PACKAGE_DIR}/LICENSE")

  if(NOT DEFINED CPACK_SYSTEM_NAME)
    set(CPACK_SYSTEM_NAME ${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR})
  endif()

  if(${CPACK_SYSTEM_NAME} MATCHES "Windows")
    if(CMAKE_CL_64)
      set(CPACK_SYSTEM_NAME win64-${CMAKE_SYSTEM_PROCESSOR})
    else()
      set(CPACK_SYSTEM_NAME win32-${CMAKE_SYSTEM_PROCESSOR})
    endif()
  endif()

  if(NOT DEFINED CPACK_PACKAGE_FILE_NAME)
    set(CPACK_PACKAGE_FILE_NAME
        "${CPACK_SOURCE_PACKAGE_FILE_NAME}-${CPACK_SYSTEM_NAME}")
  endif()

  set(CPACK_PACKAGE_EXECUTABLES "CGAL" "CGAL")

  if(WIN32 AND NOT UNIX)
    set(CPACK_GENERATOR "NSIS")
    # There is a bug in NSI that does not handle full unix paths properly. Make
    # sure there is at least one set of four (4) backlasshes.
    #set(CPACK_PACKAGE_ICON "${CGAL_INSTALLATION_PACKAGE_DIR}\\\\cgal_install.gif")
    #set(CPACK_NSIS_INSTALLED_ICON_NAME "bin\\\\CGAL.exe")
    set(CPACK_NSIS_DISPLAY_NAME
        "${CPACK_PACKAGE_INSTALL_DIRECTORY} Computational Geometry Algorithms Library"
    )
    set(CPACK_NSIS_HELP_LINK "http:\\\\\\\\www.cgal.org")
    set(CPACK_NSIS_URL_INFO_ABOUT "http:\\\\\\\\www.cgal.org")
    set(CPACK_NSIS_CONTACT "info@cgal.org")
    set(CPACK_NSIS_MODIfY_PATH ON)
  else()
    set(CPACK_STRIP_FILES "")
    set(CPACK_SOURCE_STRIP_FILES "")
  endif()

  include(CPack)

endif()

if(CGAL_BRANCH_BUILD)
  option(
    CGAL_ENABLE_CHECK_HEADERS
    "Enable the special targets \"check_pkg_headers\", \"packages_dependencies\", and \"check_pkg_<package>_headers\" and \"pkg_<package>_deps\" for each package"
    FALSE)
  unset(CGAL_COMPUTE_DEPENDENCIES CACHE) # remove that variable, no longer used

  if(CGAL_ENABLE_CHECK_HEADERS)

    option(CGAL_COPY_DEPENDENCIES
           "Copy package dependencies in source directories.
  Note that this option will modify the source directory!" FALSE)

    message("== Setting header checking ==")
    find_package(GMP REQUIRED)
    find_package(Doxygen REQUIRED)
    find_package(Eigen3 REQUIRED)
    find_package(Qt6 COMPONENTS Core Widgets OpenGL Gui REQUIRED)
    find_package(VTK COMPONENTS ImagingGeneral IOImage NO_MODULE)
    if(VTK_FOUND)
      get_target_property(VTK_INCLUDE_DIRS VTK::IOImage INTERFACE_INCLUDE_DIRECTORIES)
    endif()
    find_package(IPE)
    find_package(RS3)
    find_package(LEDA)
    find_package(OpenMesh)
    find_package(OpenCV QUIET)

    set(compile_options
        "\
${CMAKE_CXX_FLAGS} -DCGAL_EIGEN3_ENABLED -DCGAL_PROFILE \
${Qt6Widgets_DEFINITIONS} ${Qt6OpenGL_DEFINITIONS} ${Qt6Gui_DEFINITIONS} \
${Qt6OpenGL_EXECUTABLE_COMPILE_FLAGS} -fPIC \
${Qt6Gui_EXECUTABLE_COMPILE_FLAGS} \
${CGAL_3RD_PARTY_DEFINITIONS} ${CGAL_Qt6_3RD_PARTY_DEFINITIONS} \
${CGAL_DEFINITIONS}" -std=c++17)
    message("COMPILATION OPTIONS ARE :  ${compile_options}")

    if(NOT RS_FOUND AND NOT RS3_FOUND)
      set(compile_options
          "${compile_options} \
-DCGAL_ALGEBRAIC_KERNEL_RS_GMPZ_D_1=1 \
-DCGAL_ALGEBRAIC_KERNEL_RS_GMPQ_D_1=1")
      message(STATUS "Skip RS headers \"CGAL/Algebraic_kernel_rs_gmpq_d_1.h\" \
and \"CGAL/Algebraic_kernel_rs_gmpz_d_1.h\" because RS_FOUND is false.")
    else()
      set(compile_options
          "${compile_options} \
-DCGAL_USE_MPFI=1 \
-DCGAL_USE_RS=1")
    endif()
    if(NOT OpenMesh_FOUND)
      set(compile_options
          "${compile_options} \
-DCGAL_BOOST_GRAPH_GRAPH_TRAITS_POLYMESH_ARRAYKERNELT_H \
-DCGAL_BOOST_GRAPH_GRAPH_TRAITS_TRIMESH_ARRAYKERNELT_H \
-DCGAL_PROPERTIES_POLYMESH_ARRAYKERNELT_H \
-DCGAL_PROPERTIES_TRIMESH_ARRAYKERNELT_H")
      message(
        STATUS
          "Skip headers \
\"CGAL/boost/graph/graph_traits_PolyMesh_ArrayKernelT.h\", \
\"CGAL/boost/graph/graph_traits_TriMesh_ArrayKernelT.h\", \
\"CGAL/boost/graph/properties_PolyMesh_ArrayKernelT.h\", and \
\"CGAL/boost/graph/properties_TriMesh_ArrayKernelT.h\", \
because OpenMesh_FOUND is false.")
    endif()
    if(NOT OpenCV_FOUND)
      set(compile_options "${compile_options} \
-DCGAL_CLASSIFICATION_OPENCV_RANDOM_FOREST_CLASSIFIER_H")
      message(STATUS "Skip header \
\"CGAL/Classification/OpenCV_random_forest_classifier.h\" \
because OpenCV_FOUND is false")
    endif()
    if(LEDA_FOUND)
      set(compile_options "${compile_options} -DCGAL_USE_LEDA")
    endif()
    if(NOT IPE_FOUND)
      set(compile_options "${compile_options} -DCGAL_IPELET_BASE_H=1")
      message(STATUS "Skip header \"CGAL/CGAL_ipelet_base.h\" \
because IPE_FOUND is false.")
    endif()

    if(CGAL_ENABLE_CHECK_HEADERS)
      set(flag "-fsyntax-only")
    else()
      set(flag "-E")
    endif()
    if(NOT DEFINED CGAL_CHECK_SYNTAX_ONLY)
      execute_process(
        COMMAND
          ${CMAKE_CXX_COMPILER} -x c++ ${flag}
          ${CGAL_MODULES_DIR}/config/support/test_syntaxonly.cpp ERROR_QUIET
        RESULT_VARIABLE ok)
      if(ok EQUAL 0)
        set(CGAL_CHECK_SYNTAX_ONLY
            ON
            CACHE INTERNAL "")
      else()
        set(CGAL_CHECK_SYNTAX_ONLY
            OFF
            CACHE INTERNAL "")
      endif()
    endif(NOT DEFINED CGAL_CHECK_SYNTAX_ONLY)

    if(NOT CGAL_CHECK_SYNTAX_ONLY)
      message(FATAL_ERROR "Your compiler does not seem to support ${flag}.
You must disable CGAL_ENABLE_CHECK_HEADERS.")
    endif()

    ## Fill the variable include_options with all the -I and -isystem options
    set(include_options)
    foreach(incdir ${CGAL_INCLUDE_DIRS})
      list(APPEND include_options "-I${incdir}")
    endforeach()
    foreach(
      incdir
      ${VTK_INCLUDE_DIRS}
      ${LEDA_INCLUDE_DIR}
      ${OPENMESH_INCLUDE_DIR}
      ${OpenCV_INCLUDE_DIRS}
      ${RS_INCLUDE_DIR}
      ${EIGEN3_INCLUDE_DIR}
      ${GMP_INCLUDE_DIR}
      ${Qt6OpenGL_INCLUDE_DIRS}
      ${Qt6Widgets_INCLUDE_DIRS}
      ${Qt6Gui_DEFINITIONS}
      ${CGAL_3RD_PARTY_INCLUDE_DIRS}
      ${CGAL_Qt6_3RD_PARTY_INCLUDE_DIRS})
      list(APPEND include_options "-I${incdir}")
    endforeach()
    include_directories(SYSTEM ${CGAL_3RD_PARTY_INCLUDE_DIRS})

    # Build the doc
    set(DOC_DIR "${CGAL_BINARY_DIR}/build_doc")
    file(MAKE_DIRECTORY "${DOC_DIR}")
    execute_process(
      COMMAND
        "${CMAKE_COMMAND}" -DCGAL_BUILD_THREE_DOC=TRUE
        -DDOXYGEN_EXECUTABLE=${DOXYGEN_EXECUTABLE} -DBE_QUIET=TRUE
        "${CGAL_SOURCE_DIR}/Documentation/doc"
      WORKING_DIRECTORY "${DOC_DIR}")
    execute_process(
      COMMAND
        "${CMAKE_COMMAND}" --build . --target "doc"
      WORKING_DIRECTORY "${DOC_DIR}" )

  #Get the list of the documented headers
  file(GLOB html_files RELATIVE "${DOC_DIR}/doc_output/" "${DOC_DIR}/doc_output/*/*.html")
  file(GLOB example_files RELATIVE "${CGAL_SOURCE_DIR}/" "${CGAL_SOURCE_DIR}/*/examples/*/*.cpp")
  list(SORT example_files)

  find_program(AWK awk)
  set(awk_arguments [=[{ match($0, /# *include *(&lt;|[<"])(CGAL\/[^>&"]*)([>"]|&gt;)| (CGAL\/[^>&"]*\.h)/,arr); if(arr[2]!="") print arr[2]; if(arr[4]!="") print arr[4] }]=])
  message("listing headers from html files")
  foreach(f ${html_files})
    execute_process(COMMAND "${AWK}" "${awk_arguments}" "${DOC_DIR}/doc_output/${f}"
      OUTPUT_VARIABLE tmp_list)
    if (NOT "${tmp_list}" STREQUAL "")
      string(REPLACE "\n" ";" tmp_list ${tmp_list})
      LIST( APPEND list_of_documented_headers ${tmp_list})
    endif()
  endforeach()

  message("listing headers from examples files")
  foreach(f ${example_files})
    execute_process(COMMAND "${AWK}" "${awk_arguments}" "${CGAL_SOURCE_DIR}/${f}"
      OUTPUT_VARIABLE tmp_list)
    if (NOT "${tmp_list}" STREQUAL "")
      string(REPLACE "\n" ";" tmp_list ${tmp_list})
      LIST( APPEND list_of_documented_headers ${tmp_list})
    endif()
  endforeach()
  message("removing duplicates:")
  list(REMOVE_DUPLICATES list_of_documented_headers)
  message("sorting:")
  list(SORT list_of_documented_headers)
  #string(REPLACE ";" " \n" list_of_documented_headers "${list_of_documented_headers}")
  list(SORT list_of_documented_headers)
  list(REMOVE_DUPLICATES list_of_documented_headers)

  include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules/list_of_whitelisted_headers.cmake
      OPTIONAL RESULT_VARIABLE has_list_of_whitelisted_headers)
    message("list_of_whitelisted_headers: ${list_of_whitelisted_headers}")
    ## Loop on package and headers
    set(check_pkg_target_list)
    #get build dir for removal from deps
    get_filename_component(BUILD_DIR ${CGAL_BINARY_DIR} NAME)
    foreach(package ${CGAL_CONFIGURED_PACKAGES_NAMES})
      if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/../${package}/include)
        set(check_pkg_headers_depends "")
        file(GLOB_RECURSE ${package}_HEADERS #don't use a hard coded list because the ones for the documented headers is enough.
          RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}/../${package}/include"
          "${CMAKE_CURRENT_SOURCE_DIR}/../${package}/include/CGAL/*.h*")

        list(SORT ${package}_HEADERS)
        foreach(header ${${package}_HEADERS})
          set(skip_hdr FALSE)
          if(NOT header IN_LIST list_of_documented_headers)
#            message(STATUS "Skip non-documented header \"${header}\".")
            set(skip_hdr TRUE)
          endif()
          if(NOT VTK_FOUND)
            string(REGEX MATCH ".*vtk.*" is_a_vtk_header ${header})
            if(is_a_vtk_header)
              message(STATUS "Skip VTK header \"${header}\" because \
VTK_FOUND is false.")
              set(skip_hdr TRUE)
            endif()
          endif()
          string(REGEX MATCH "CGAL/leda_.*" is_a_leda_header ${header})
          if(NOT LEDA_FOUND AND is_a_leda_header)
            message(STATUS "Skip LEDA header \"${header}\" because \
LEDA_FOUND is false.")
            set(skip_hdr TRUE)
          endif()
          if(has_list_of_whitelisted_headers
              AND header IN_LIST list_of_whitelisted_headers)
            message(STATUS "Skip whitelisted header \"${header}\".")
            set(skip_hdr TRUE)
          endif()
          if(NOT skip_hdr)
            list(APPEND list_of_headers_to_test ${header})
            string(REPLACE "/" "__" header2 "${header}")
            string(REPLACE "." "_" header2 "${header2}")
            string(REPLACE ";" " " include_options_str "${include_options}")
            string(REPLACE ";" " " compile_options_str "${compile_options}")
            separate_arguments(
              CMD UNIX_COMMAND "${CMAKE_CXX_COMPILER} ${compile_options_str}
${include_options_str} -x c++ ${flag} -H \
${CMAKE_CURRENT_SOURCE_DIR}/../${package}/include/${header}"
              # The header Algebraic_kernel_rs_gmpz_d_1.h is skipped on purpose: it
              # depends on RS.
            )
            # CMD2 is CMD without the -H option
            separate_arguments(
              CMD2 UNIX_COMMAND "${CMAKE_CXX_COMPILER} ${compile_options_str}
${include_options_str} -x c++ ${flag} \
${CMAKE_CURRENT_SOURCE_DIR}/../${package}/include/${header}"
              # The header Algebraic_kernel_rs_gmpz_d_1.h is skipped on purpose: it
              # depends on RS.
            )
            set(chk_header_name
                ${CGAL_BINARY_DIR}/package_info/${package}/check_headers/check_${header2}
            )
            add_custom_command(
              OUTPUT ${chk_header_name}
              COMMAND
                ${CMAKE_COMMAND} -DCERR:STRING=${chk_header_name}
                "-DCMD:STRING=${CMD}" "-DCMD2:STRING=${CMD2}" -P
                "${CGAL_MODULES_DIR}/run_cmd_redirection_cerr.cmake"
              DEPENDS
                ${CMAKE_CURRENT_SOURCE_DIR}/../${package}/include/${header}
              DEPENDS ${CGAL_MODULES_DIR}/run_cmd_redirection_cerr.cmake
              VERBATIM
              COMMENT "Check header ${package}/include/${header}")
            list(APPEND check_pkg_headers_depends ${chk_header_name})
            if(${header2}_deps)
              list(REMOVE_DUPLICATES ${header2}_deps)
            endif()
          endif() # end the if on list_of_documented_headers
        endforeach() # look on headers
        add_custom_target(check_pkg_${package}_headers
                          DEPENDS ${check_pkg_headers_depends})
        add_custom_command(
          OUTPUT ${CGAL_BINARY_DIR}/package_info/${package}/included_headers
                 ${CGAL_BINARY_DIR}/package_info/${package}/dependencies
          DEPENDS ${check_pkg_headers_depends}
          COMMENT "Compute dependencies of ${package}"
          COMMAND
            ${CMAKE_COMMAND} -DCGAL_PACKAGES_PREFIX=${CGAL_SOURCE_DIR}
            -DCGAL_COPY_DEPENDENCIES=${CGAL_COPY_DEPENDENCIES}
            -DCGAL_COPY_PATH=${CGAL_SOURCE_DIR}/${package}/package_info/${package}/dependencies
            -DOUTPUT_HEADERS_LIST=${CGAL_BINARY_DIR}/package_info/${package}/included_headers
            -DOUTPUT_PACKAGES_LIST=${CGAL_BINARY_DIR}/package_info/${package}/dependencies
            -P "${CGAL_MODULES_DIR}/process_dependencies.cmake"
            ${check_pkg_headers_depends})
        add_custom_target(
          pkg_${package}_deps
          DEPENDS ${CGAL_BINARY_DIR}/package_info/${package}/dependencies)
        list(APPEND check_pkg_target_list check_pkg_${package}_headers)
        list(APPEND packages_deps
             ${CGAL_BINARY_DIR}/package_info/${package}/dependencies)
      endif() # if the package has an include directory
      if(packages_deps)
        list(REMOVE_DUPLICATES packages_deps)
      endif()
    endforeach() # loop on packages
    #Now check that a cpp file including all documented headers compiles
    file(WRITE ${CGAL_BINARY_DIR}/test_headers.cpp "#include <CGAL/Periodic_3_mesh_3/config.h>\n")
    foreach(header ${list_of_headers_to_test})
      file(APPEND ${CGAL_BINARY_DIR}/test_headers.cpp "#include <${header}>\n")
    endforeach() #loop on headers to include in test file
    file(APPEND ${CGAL_BINARY_DIR}/test_headers.cpp "int main(){}\n")

    add_custom_target(check_headers DEPENDS ${check_pkg_target_list})

    string(REPLACE " " ";" compile_options ${compile_options})
    string(REPLACE " " ";" include_options ${include_options_str})
    add_custom_command(
      TARGET check_headers
      POST_BUILD
      COMMAND ${CMAKE_CXX_COMPILER} ${compile_options} ${include_options} -x
              c++ -fsyntax-only ${CGAL_BINARY_DIR}/test_headers.cpp)

    add_custom_target(packages_dependencies DEPENDS ${packages_deps})
    message(
      "  \n\
  You can now check the headers with the target `check_headers`\n\
  and the package dependencies with the target `packages_dependencies`.\n\
  Results are in the `package_info/` sub-directory of the build directory:\n\
  \n\
    - package_info/<package>/dependencies\n\
    - package_info/<package>/included_headers\n\
    - package_info/<package>/check_headers/ (error messages from \
the headers checks)\n")
    message("== Setting header checking (DONE) ==\n")
  endif()
endif(CGAL_BRANCH_BUILD)

if(NOT CGAL_BRANCH_BUILD AND EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/doc")
  # in a non-branch build this is the top-level CMakeLists.txt
  add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/doc")
endif()

#print some info about versions
if(RUNNING_CGAL_AUTO_TEST OR CGAL_TEST_SUITE)
  find_package(Qt6 QUIET COMPONENTS Core)
  find_package(Boost)
  if(NOT CGAL_DISABLE_GMP)
    find_package(GMP)
    find_package(MPFR)
    get_dependency_version(GMP)
    get_dependency_version(MPFR)
  elseif(WITH_LEDA)#CGAL_DISABLE_GMP
    find_package(LEDA)
  endif()#NOT CGAL_DISABLE_GMP
  message(
    STATUS
    "USING BOOST_VERSION = '${Boost_MAJOR_VERSION}.${Boost_MINOR_VERSION}.${Boost_SUBMINOR_VERSION}'"
    )
  if(Qt6_FOUND)
    message(STATUS "USING Qt6_VERSION = '${Qt6Core_VERSION}'")
  endif()#Qt6_FOUND
endif()#RUNNING_CGAL_AUTO_TEST OR CGAL_TEST_SUITE
