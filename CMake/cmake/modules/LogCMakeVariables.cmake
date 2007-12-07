# Generic and CGAL-specific CMake variable logging
# (see http://www.cmake.org/Wiki/CMake_Useful_Variables/Logging_Useful_Variables)


# Utility macro: log a (text) variable
MACRO(LOG_VARIABLE VAR_NAME)
    MESSAGE(STATUS "${VAR_NAME}: " ${${VAR_NAME}})
ENDMACRO(LOG_VARIABLE)


# Utility macro: log a list
MACRO(LOG_LIST LIST_NAME SEPARATOR)
    SET(_expanded_LIST)
    FOREACH(_current ${${LIST_NAME}})
        SET(_expanded_LIST "${_expanded_LIST} ${SEPARATOR}${_current}")
    ENDFOREACH(_current ${${LIST_NAME}})
    MESSAGE(STATUS "${LIST_NAME}: " ${_expanded_LIST})
ENDMACRO(LOG_LIST)


# Log all generic and CGAL-specific CMake variable logging
MACRO (LOG_CMAKE_VARIABLES COMMENT)

    MESSAGE(STATUS "--------------- Start CMake Variable Logging (${COMMENT}) ---------------")

    # this is the file name and line number of the file where this variable is used.
    LOG_VARIABLE(CMAKE_CURRENT_LIST_FILE)
    LOG_VARIABLE(CMAKE_CURRENT_LIST_LINE)


    # if you are building in-source, this is the same as CMAKE_SOURCE_DIR, otherwise
    # this is the top level directory of your build tree
    MESSAGE(STATUS "CMAKE_BINARY_DIR:         " ${CMAKE_BINARY_DIR})
    MESSAGE(STATUS "CGAL_BINARY_DIR:          " ${CGAL_BINARY_DIR})

    # if you are building in-source, this is the same as CMAKE_CURRENT_SOURCE_DIR, otherwise this
    # is the directory where the compiled or generated files from the current CMakeLists.txt will go to
    MESSAGE(STATUS "CMAKE_CURRENT_BINARY_DIR: " ${CMAKE_CURRENT_BINARY_DIR})

    # this is the directory, from which cmake was started, i.e. the top level source directory
    MESSAGE(STATUS "CMAKE_SOURCE_DIR:         " ${CMAKE_SOURCE_DIR})
    MESSAGE(STATUS "CGAL_SOURCE_DIR:          " ${CGAL_SOURCE_DIR})

    # this is the directory where the currently processed CMakeLists.txt is located in
    MESSAGE(STATUS "CMAKE_CURRENT_SOURCE_DIR: " ${CMAKE_CURRENT_SOURCE_DIR})

    # contains the full path to the top level directory of your build tree
    MESSAGE(STATUS "PROJECT_BINARY_DIR:       " ${PROJECT_BINARY_DIR})

    # contains the full path to the root of your project source directory,
    # i.e. to the nearest directory where CMakeLists.txt contains the PROJECT() command
    MESSAGE(STATUS "PROJECT_SOURCE_DIR:       " ${PROJECT_SOURCE_DIR})

    # set this variable to specify a common place where CMake should put all executable files
    # (instead of CMAKE_CURRENT_BINARY_DIR)
    MESSAGE(STATUS "EXECUTABLE_OUTPUT_PATH:   " ${EXECUTABLE_OUTPUT_PATH})

    # set this variable to specify a common place where CMake should put all libraries
    # (instead of CMAKE_CURRENT_BINARY_DIR)
    MESSAGE(STATUS "LIBRARY_OUTPUT_PATH:      " ${LIBRARY_OUTPUT_PATH})

    # tell CMake to search first in directories listed in CMAKE_MODULE_PATH
    # when you use FIND_PACKAGE() or INCLUDE()
    LOG_LIST(CMAKE_MODULE_PATH " ")


    # Configure install locations.
    MESSAGE(STATUS "CMAKE_INSTALL_PREFIX      " ${CMAKE_INSTALL_PREFIX})
    MESSAGE(STATUS "CGAL_INSTALL_DIR          " ${CGAL_INSTALL_DIR})
    MESSAGE(STATUS "CGAL_LIB_INSTALL_DIR      " ${CGAL_LIB_INSTALL_DIR})
    MESSAGE(STATUS "CGAL_INCLUDE_INSTALL_DIR  " ${CGAL_INCLUDE_INSTALL_DIR})
    MESSAGE(STATUS "CGAL_BIN_INSTALL_DIR      " ${CGAL_BIN_INSTALL_DIR})
    MESSAGE(STATUS "CGAL_AUXILIARY_INSTALL_DIR" ${CGAL_AUXILIARY_INSTALL_DIR})
    MESSAGE(STATUS "CGAL_MODULE_INSTALL_DIR   " ${CGAL_MODULE_INSTALL_DIR})


    # this is the complete path of the cmake which runs currently (e.g. /usr/local/bin/cmake)
    LOG_VARIABLE(CMAKE_COMMAND)

    # this is CMake's installation directory
    LOG_VARIABLE(CMAKE_ROOT)

    # this is used when searching for include files e.g. using the FIND_PATH() command.
    LOG_VARIABLE(CMAKE_INCLUDE_PATH)

    # this is used when searching for libraries e.g. using the FIND_LIBRARY() command.
    LOG_VARIABLE(CMAKE_LIBRARY_PATH)

    # the complete system name, e.g. "Linux-2.4.22", "FreeBSD-5.4-RELEASE" or "Windows 5.1"
    LOG_VARIABLE(CMAKE_SYSTEM)

    # the short system name, e.g. "Linux", "FreeBSD" or "Windows"
    LOG_VARIABLE(CMAKE_SYSTEM_NAME)

    # only the version part of CMAKE_SYSTEM
    LOG_VARIABLE(CMAKE_SYSTEM_VERSION)

    # the processor name (e.g. "Intel(R) Pentium(R) M processor 2.00GHz")
    LOG_VARIABLE(CMAKE_SYSTEM_PROCESSOR)

    # is TRUE on all UNIX-like OS's, including Apple OS X and CygWin
    LOG_VARIABLE(UNIX)

    # is TRUE on Windows, including CygWin
    LOG_VARIABLE(WIN32)

    # is TRUE on Apple OS X
    LOG_VARIABLE(APPLE)

    # is TRUE when using the MinGW compiler in Windows
    LOG_VARIABLE(MINGW)

    # is TRUE on Windows when using the CygWin version of cmake
    LOG_VARIABLE(CYGWIN)

    # is TRUE on Windows when using a Borland compiler
    LOG_VARIABLE(BORLAND)

    # Microsoft compiler
    LOG_VARIABLE(MSVC)
    LOG_VARIABLE(MSVC_IDE)
    LOG_VARIABLE(MSVC60)
    LOG_VARIABLE(MSVC70)
    LOG_VARIABLE(MSVC71)
    LOG_VARIABLE(MSVC80)
    LOG_VARIABLE(CMAKE_COMPILER_2005)


    # set this to true if you don't want to rebuild the object files if the rules have changed,
    # but not the actual source files or headers (e.g. if you changed the some compiler switches)
    LOG_VARIABLE(CMAKE_SKIP_RULE_DEPENDENCY)

    # since CMake 2.1 the install rule depends on all, i.e. everything will be built before installing.
    # If you don't like this, set this one to true.
    LOG_VARIABLE(CMAKE_SKIP_INSTALL_ALL_DEPENDENCY)

    # If set, runtime paths are not added when using shared libraries. Default it is set to OFF
    LOG_VARIABLE(CMAKE_SKIP_RPATH)

    # set this to true if you are using makefiles and want to see the full compile and link
    # commands instead of only the shortened ones
    LOG_VARIABLE(CMAKE_VERBOSE_MAKEFILE)

    # this will cause CMake to not put in the rules that re-run CMake. This might be useful if
    # you want to use the generated build files on another machine.
    LOG_VARIABLE(CMAKE_SUPPRESS_REGENERATION)


    # A simple way to get switches to the compiler is to use ADD_DEFINITIONS().
    # But there are also two variables exactly for this purpose:

    # the compiler flags for compiling C sources
#    LOG_VARIABLE(CMAKE_C_FLAGS)

    # the compiler flags for compiling C++ sources
    LOG_VARIABLE(CMAKE_CXX_FLAGS)

    # Include directories
    GET_DIRECTORY_PROPERTY(INCLUDE_DIRECTORIES INCLUDE_DIRECTORIES)
    LOG_LIST(INCLUDE_DIRECTORIES "-I")

    # Choose the type of build.  Example: SET(CMAKE_BUILD_TYPE Debug)
    LOG_VARIABLE(CMAKE_BUILD_TYPE)

    # if this is set to ON, then all libraries are built as shared libraries by default.
    LOG_VARIABLE(BUILD_SHARED_LIBS)

    # the compiler used for C files
#    LOG_VARIABLE(CMAKE_C_COMPILER)

    # the compiler used for C++ files
    LOG_VARIABLE(CMAKE_CXX_COMPILER)

    # if the compiler is a variant of gcc, this should be set to 1
#     LOG_VARIABLE(CMAKE_COMPILER_IS_GNUCC)

    # if the compiler is a variant of g++, this should be set to 1
    LOG_VARIABLE(CMAKE_COMPILER_IS_GNUCXX)

    # the tools for creating libraries
    LOG_VARIABLE(CMAKE_AR)
    LOG_VARIABLE(CMAKE_RANLIB)

#     # Rules for C++ sources:
#     LOG_VARIABLE(CMAKE_CXX_CREATE_SHARED_LIBRARY)
#     LOG_VARIABLE(CMAKE_CXX_CREATE_SHARED_MODULE)
#     LOG_VARIABLE(CMAKE_CXX_CREATE_STATIC_LIBRARY)
#     LOG_VARIABLE(CMAKE_CXX_COMPILE_OBJECT)
#     LOG_VARIABLE(CMAKE_CXX_LINK_EXECUTABLE)

    # CGAL Platform configuration tests
    LOG_VARIABLE(CGAL_CFG_BOOL_IN_TEMPLATE_BUG)
    LOG_VARIABLE(CGAL_CFG_CCTYPE_MACRO_BUG)
    LOG_VARIABLE(CGAL_CFG_COMMA_BUG)
    LOG_VARIABLE(CGAL_CFG_CONVERSION_OPERATOR_BUG)
    LOG_VARIABLE(CGAL_CFG_DEEP_DEPENDENT_TEMPLATE_BUG)
    LOG_VARIABLE(CGAL_CFG_DENORMALS_COMPILE_BUG)
    LOG_VARIABLE(CGAL_CFG_IEEE_754_BUG)
    LOG_VARIABLE(CGAL_CFG_ISTREAM_INT_BUG)
    LOG_VARIABLE(CGAL_CFG_LONGNAME_BUG)
    LOG_VARIABLE(CGAL_CFG_MATCHING_BUG_3)
    LOG_VARIABLE(CGAL_CFG_MATCHING_BUG_4)
    LOG_VARIABLE(CGAL_CFG_MATCHING_BUG_5)
    LOG_VARIABLE(CGAL_CFG_MATCHING_BUG_6)
    LOG_VARIABLE(CGAL_CFG_MISSING_TEMPLATE_VECTOR_CONSTRUCTORS_BUG)
    LOG_VARIABLE(CGAL_CFG_NESTED_CLASS_FRIEND_DECLARATION_BUG)
    LOG_VARIABLE(CGAL_CFG_NET2003_MATCHING_BUG)
    LOG_VARIABLE(CGAL_CFG_NO_KOENIG_LOOKUP)
    LOG_VARIABLE(CGAL_CFG_NO_LIMITS)
    LOG_VARIABLE(CGAL_CFG_NO_LOCALE)
    LOG_VARIABLE(CGAL_CFG_NO_LONG_DOUBLE_IO)
    LOG_VARIABLE(CGAL_CFG_NO_LONG_LONG)
    LOG_VARIABLE(CGAL_CFG_NO_NEXTAFTER)
    LOG_VARIABLE(CGAL_CFG_NO_STDC_NAMESPACE)
    LOG_VARIABLE(CGAL_CFG_NO_STL)
    LOG_VARIABLE(CGAL_CFG_NO_TMPL_IN_TMPL_DEPENDING_FUNCTION_PARAM)
    LOG_VARIABLE(CGAL_CFG_NO_TMPL_IN_TMPL_PARAM)
    LOG_VARIABLE(CGAL_CFG_NO_TWO_STAGE_NAME_LOOKUP)
    LOG_VARIABLE(CGAL_CFG_NUMERIC_LIMITS_BUG)
    LOG_VARIABLE(CGAL_CFG_OUTOFLINE_MEMBER_DEFINITION_BUG)
    LOG_VARIABLE(CGAL_CFG_OUTOFLINE_TEMPLATE_MEMBER_DEFINITION_BUG_2)
    LOG_VARIABLE(CGAL_CFG_OUTOFLINE_TEMPLATE_MEMBER_DEFINITION_BUG)
    LOG_VARIABLE(CGAL_CFG_SUNPRO_RWSTD)
    LOG_VARIABLE(CGAL_CFG_TYPENAME_BEFORE_DEFAULT_ARGUMENT_BUG)
    LOG_VARIABLE(CGAL_CFG_USING_BASE_MEMBER_BUG_2)
    LOG_VARIABLE(CGAL_CFG_USING_BASE_MEMBER_BUG_3)
    LOG_VARIABLE(CGAL_CFG_USING_BASE_MEMBER_BUG)


    # CGAL third party libraries
    LOG_VARIABLE(CGAL_USE_BOOST)
    LOG_VARIABLE(CGAL_USE_BOOST_PROGRAM_OPTIONS)
    LOG_VARIABLE(CGAL_USE_X11)
    LOG_VARIABLE(CGAL_USE_GMP)
    LOG_VARIABLE(CGAL_USE_GMPXX)
    LOG_VARIABLE(CGAL_USE_MPFR)
    LOG_VARIABLE(CGAL_USE_CORE)
    LOG_VARIABLE(CGAL_USE_CGAL_CORE)
    LOG_VARIABLE(CGAL_USE_ZLIB)
    LOG_VARIABLE(CGAL_USE_LIDIA)
    LOG_VARIABLE(CGAL_USE_LEDA)
    LOG_VARIABLE(CGAL_USE_LEDAWIN)
    LOG_VARIABLE(CGAL_USE_QT)
    LOG_VARIABLE(CGAL_USE_TAUCS)
    LOG_LIST(CGAL_3RD_PARTY_DEFINITIONS " ")
    LOG_LIST(CGAL_3RD_PARTY_INCLUDE_DIRS "-I")
    LOG_LIST(CGAL_3RD_PARTY_LIBRARIES " ")


    # Settings to compile CGAL libraries and with CGAL libraries
    LOG_LIST(CGAL_DEFINITIONS " ")
    LOG_LIST(CGAL_INCLUDE_DIRS "-I")
    LOG_LIST(CGAL_LIB_CXX_FLAGS " ")
    LOG_VARIABLE(CGAL_LIB_EXE_LINKER_FLAGS)
    LOG_VARIABLE(CGAL_LIB_SHARED_LINKER_FLAGS)
    LOG_LIST(CGAL_LIBRARIES " ")


    MESSAGE(STATUS "--------------- End CMake Variable Logging (${COMMENT}) ---------------")

ENDMACRO (LOG_CMAKE_VARIABLES)
