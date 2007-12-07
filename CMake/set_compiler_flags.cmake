# set_compiler_flags.cmake
# Compute special compiler flags required by CGAL on some platforms.
#
# The code below is translated from install_cgal's set_compiler_flags()
# and uses the same variable names.
#
# TODO: On each platform: Test and adapt this code. Some settings are already provided by CMake.
# Done: - Linux/gcc 3.4
#       - Windows/Visual C++ 7.1
# Note: Some of the compilers here are probably not supported by CMake.


#
# Default value for the compiler/os specific variables set here:
#

# Common compiler settings
SET(ADDITIONAL_CXXFLAGS "")

# Common linker settings
SET(ADDITIONAL_LDFLAGS "")

# Extra settings when compiling CGAL libs (static)
SET(CGAL_LIB_CXXFLAGS "")

# Extra settings when linking CGAL libs (static)
SET(CGAL_LIB_LDFLAGS "")

# Linker name (for static libs).
# Obsolete: built-in in CMake.
# SET(CGAL_LIB_CREATE "ar cr ''")

# Extra settings when compiling CGAL libs (dynamic).
# Obsolete: -fpic option is built-in in CMake.
#SET(CGAL_SHARED_LIB_CXXFLAGS "-fpic")

# Extra settings when linking CGAL libs (dynamic)
SET(CGAL_SHARED_LIB_LDFLAGS "")

# Linker name (for dynamic libs).
# Obsolete: built-in in CMake.
# SET(CGAL_SHARED_LIB_CREATE "${CMAKE_CXX_COMPILER} -shared"

# Rpath setting when linking CGAL libs (dynamic).
# Obsolete: CMake provides rpath control through CMAKE_SKIP_RPATH.
# SET(RUNTIME_LINKER_FLAG '')

# Math lib when linking CGAL libs (static).
# Seems built-in in CMake.
# SET(MATH_LIB "m")


#
# Set values based on current compiler/os specific:
#

# Get compiler's version
EXEC_PROGRAM(${CMAKE_CXX_COMPILER} ARGS --version
             OUTPUT_VARIABLE CMAKE_CXX_COMPILER_VERSION)

if (BORLAND)
    #### settings Borland C++ compiler
    SET(ADDITIONAL_CXXFLAGS "-P -vi- -w-inl")
    SET(CGAL_LIB_CXXFLAGS "-P")
#     SET(CGAL_LIB_CREATE "tlib /C /P512 ''")
#     SET(MATH_LIB "")
endif (BORLAND)

if (CMAKE_CXX_COMPILER MATCHES "^(ICL|icl)")
    #### settings Intel C++ compiler
    SET(ADDITIONAL_CXXFLAGS "-Op -TP -GR -EHsc -Zm900 -nologo -Zc:wchar_t -Zc:forScope -MD")
    SET(ADDITIONAL_LDFLAGS "-nologo")
#     SET(CGAL_LIB_CREATE "LIB /OUT:")
    SET(CGAL_SHARED_LIB_CXXFLAGS "")
#     SET(CGAL_SHARED_LIB_CREATE "")
#     SET(MATH_LIB "")
endif (CMAKE_CXX_COMPILER MATCHES "^(ICL|icl)")

if (CMAKE_CXX_COMPILER MATCHES "^(CL|cl)")
    #### settings Microsoft Visual C++ compiler
#    SET(ADDITIONAL_CXXFLAGS "-TP -GR -EHsc -Zm900 -nologo -MD")
    SET(ADDITIONAL_CXXFLAGS "-TP -GR -EHsc -nologo -MD") # CMake sets -Zm1000
    SET(ADDITIONAL_LDFLAGS "-nologo")
#   SET(CGAL_LIB_CREATE "LIB /OUT:")
    SET(CGAL_SHARED_LIB_CXXFLAGS "")
#   SET(CGAL_SHARED_LIB_CREATE "")
#   SET(MATH_LIB "")
endif (CMAKE_CXX_COMPILER MATCHES "^(CL|cl)")

if(CMAKE_SYSTEM MATCHES "IRIX.*5" AND CMAKE_CXX_COMPILER MATCHES "^(CC|cc)")
    #### settings for sgi mipspro compiler on irix5
    SET(CGAL_SHARED_LIB_CXXFLAGS "")
#     SET(CGAL_SHARED_LIB_LDFLAGS "-lm") # built-in in CMake
#     SET(RUNTIME_LINKER_FLAG "-rpath ")
endif(CMAKE_SYSTEM MATCHES "IRIX.*5" AND CMAKE_CXX_COMPILER MATCHES "^(CC|cc)")

if(CMAKE_SYSTEM MATCHES "IRIX.*6" AND
   CMAKE_CXX_COMPILER MATCHES "^(CC|cc)" AND CMAKE_CXX_COMPILER_VERSION MATCHES "7\\.3")
    #### settings for sgi mipspro compiler V7.3 on irix6
    SET(ADDITIONAL_CXXFLAGS "-LANG:std")
    SET(ADDITIONAL_LDFLAGS "-LANG:std")
#     SET(CGAL_LIB_CREATE "${CMAKE_CXX_COMPILER} -ar -o''")
    SET(CGAL_SHARED_LIB_CXXFLAGS "")
#     SET(RUNTIME_LINKER_FLAG "-rpath ")

else(CMAKE_SYSTEM MATCHES "IRIX.*6" AND
     CMAKE_CXX_COMPILER MATCHES "^(CC|cc)" AND CMAKE_CXX_COMPILER_VERSION MATCHES "7\\.3")

    if(CMAKE_SYSTEM MATCHES "IRIX.*6" AND CMAKE_CXX_COMPILER MATCHES "^(CC|cc)")
        #### settings for sgi mipspro compiler (except V7.3) on irix6
#         SET(CGAL_LIB_CREATE "${CMAKE_CXX_COMPILER} -ar -o''")
        SET(CGAL_SHARED_LIB_CXXFLAGS "")
#         SET(RUNTIME_LINKER_FLAG "-rpath ")
    endif(CMAKE_SYSTEM MATCHES "IRIX.*6" AND CMAKE_CXX_COMPILER MATCHES "^(CC|cc)")

endif(CMAKE_SYSTEM MATCHES "IRIX.*6" AND
      CMAKE_CXX_COMPILER MATCHES "^(CC|cc)" AND CMAKE_CXX_COMPILER_VERSION MATCHES "7\\.3")

if(CMAKE_SYSTEM MATCHES "IRIX" AND CMAKE_COMPILER_IS_GNUCXX)
    #### settings for gcc on irix
    SET(ADDITIONAL_CXXFLAGS "-Wall")
#     SET(CGAL_SHARED_LIB_LDFLAGS "-lm") # built-in in CMake
#     SET(RUNTIME_LINKER_FLAG "-Xlinker -rpath -Xlinker ")
endif(CMAKE_SYSTEM MATCHES "IRIX" AND CMAKE_COMPILER_IS_GNUCXX)

if(CMAKE_SYSTEM MATCHES "SunOS-5" AND CMAKE_CXX_COMPILER MATCHES "^(CC|cc)")
    #### settings for sunpro compiler on solaris
    SET(ADDITIONAL_CXXFLAGS "-features=extensions -D_RWSTD_ALLOCATOR")
#     SET(CGAL_LIB_CREATE "${CMAKE_CXX_COMPILER} -xar -o ''")
    SET(CGAL_SHARED_LIB_CXXFLAGS "-pic") # probably built-in in CMake
#     SET(CGAL_SHARED_LIB_CREATE "${CMAKE_CXX_COMPILER} -G")
#     SET(RUNTIME_LINKER_FLAG "-R ")
endif(CMAKE_SYSTEM MATCHES "SunOS-5" AND CMAKE_CXX_COMPILER MATCHES "^(CC|cc)")

if(CMAKE_SYSTEM MATCHES "SunOS-5" AND CMAKE_COMPILER_IS_GNUCXX)
    #### settings for gcc on solaris
#     SET(CGAL_SHARED_LIB_CREATE "${CMAKE_CXX_COMPILER} -G")
#     SET(RUNTIME_LINKER_FLAG "-R ")
endif(CMAKE_SYSTEM MATCHES "SunOS-5" AND CMAKE_COMPILER_IS_GNUCXX)

if(CMAKE_SYSTEM MATCHES "alpha.*Linux" AND CMAKE_COMPILER_IS_GNUCXX)
    #### settings for g++ on alpha-linux (special FPU handling)
    #### LONG_NAME_PROBLEM is cured by disabling debugging
    SET(ADDITIONAL_CXXFLAGS "-Wall -mieee -mfp-rounding-mode=d")
#     SET(RUNTIME_LINKER_FLAG '-Wl,-R')
endif(CMAKE_SYSTEM MATCHES "alpha.*Linux" AND CMAKE_COMPILER_IS_GNUCXX)

if(CMAKE_SYSTEM MATCHES "^Linux" AND CMAKE_COMPILER_IS_GNUCXX)
    #### settings for gcc on linux (except alpha-linux)
    SET(ADDITIONAL_CXXFLAGS "-Wall")
#     SET(RUNTIME_LINKER_FLAG '-Wl,-R')
endif(CMAKE_SYSTEM MATCHES "^Linux" AND CMAKE_COMPILER_IS_GNUCXX)

if(CMAKE_SYSTEM MATCHES "Linux" AND CMAKE_CXX_COMPILER MATCHES "^(icc|icpc)")
    #### settings for icc on linux
    # -mp is required for correct enough floating point operations
    # necessary for interval arithmetic.
    SET(ADDITIONAL_CXXFLAGS "-mp")
#     SET(RUNTIME_LINKER_FLAG '-Wl,-R')
endif(CMAKE_SYSTEM MATCHES "Linux" AND CMAKE_CXX_COMPILER MATCHES "^(icc|icpc)")

if(CYGWIN AND CMAKE_COMPILER_IS_GNUCXX)
    #### settings for gcc on Cygwin
    SET(ADDITIONAL_CXXFLAGS "-Wall")
endif(CYGWIN AND CMAKE_COMPILER_IS_GNUCXX)

if(MINGW AND CMAKE_COMPILER_IS_GNUCXX)
    #### settings for gcc on mingw
    SET(ADDITIONAL_CXXFLAGS "-Wall")
endif(MINGW AND CMAKE_COMPILER_IS_GNUCXX)

if(APPLE AND CMAKE_COMPILER_IS_GNUCXX)
    #### settings for gcc on Darwin (MacOSX)
    SET(ADDITIONAL_CXXFLAGS "-Wall")
    SET(CGAL_SHARED_LIB_CXXFLAGS "-fno-common")
#     SET(CGAL_SHARED_LIB_CREATE "${CMAKE_CXX_COMPILER} -dynamiclib")
endif(APPLE AND CMAKE_COMPILER_IS_GNUCXX)

if (CMAKE_CXX_COMPILER MATCHES "^(mwcc|MWCC)")
    #### settings for Metrowerks C++ compiler
    SET(ADDITIONAL_CXXFLAGS "-gccincludes -lang c++ -msgstyle gcc ")
    SET(CGAL_LIB_CXXFLAGS "")
#     SET(CGAL_LIB_CREATE "mwld -library -o ")
#     SET(MATH_LIB "")
endif (CMAKE_CXX_COMPILER MATCHES "^(mwcc|MWCC)")

if(CMAKE_SYSTEM MATCHES "Linux" AND CMAKE_CXX_COMPILER MATCHES "^pgCC")
    #### settings for Portland Group Compiler on linux
#     SET(RUNTIME_LINKER_FLAG '-Wl,-R')
    # PGCC has long name problems with "-g".
endif(CMAKE_SYSTEM MATCHES "Linux" AND CMAKE_CXX_COMPILER MATCHES "^pgCC")

#### else: unknown compiler
# TODO: prompt an error if unknown compiler

# Version specific adjustments for gcc
if(CMAKE_COMPILER_IS_GNUCXX)
    # If gcc 2.x
    IF(CMAKE_CXX_COMPILER_VERSION MATCHES "2\\.[0-9]")
        SET(ADDITIONAL_CXXFLAGS "${ADDITIONAL_CXXFLAGS} -ftemplate-depth-50")
    ENDIF(CMAKE_CXX_COMPILER_VERSION MATCHES "2\\.[0-9]")

    # If gcc 3.x
    IF(CMAKE_CXX_COMPILER_VERSION MATCHES "3\\.[0-9]")
        # Nothing to add
    ELSE(CMAKE_CXX_COMPILER_VERSION MATCHES "3\\.[0-9]")
        # g++-4.x and after: we need -frounding-math to avoid FPU rounding mode optimisations
        SET(ADDITIONAL_CXXFLAGS "${ADDITIONAL_CXXFLAGS} -ftemplate-depth-50")
    ENDIF(CMAKE_CXX_COMPILER_VERSION MATCHES "3\\.[0-9]")
endif(CMAKE_COMPILER_IS_GNUCXX)

