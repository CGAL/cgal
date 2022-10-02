# - Try to find FileGDBAPI (FileGDBAPI.lib on Windows and v8_base.x64.a on Linux)
# Once done this will define
#  FileGDBAPI_FOUND - System has FileGDBAPI
#  FileGDBAPI_INCLUDE_DIR - The FileGDBAPI include directories
#  FileGDBAPI_LIBRARY - The library needed to use FileGDBAPI
#  FileGDBAPI_LIBRARY_DIR - The directory where lib files are.

set(FileGDBAPI_NAMES_RELEASE FileGDBAPI)
set(FileGDBAPI_NAMES_DEBUG FileGDBAPId ${FileGDBAPI_NAMES_RELEASE})

find_path(FileGDBAPI_INCLUDE_DIR NAMES FileGDBAPI.h
  PATH_SUFFIXES include FileGDBAPI include/FileGDBAPI
  PATHS /opt/libvFileGDBAPI-$ENV{FileGDBAPI_VER}
  HINTS ENV FileGDBAPI_INC_DIR ENV FileGDBAPI_DIR)

# CMake>=2.6 supports the notation "debug XXd optimized XX"
set(FileGDBAPI_HINTS ENV FileGDBAPI_LIB_DIR ENV FileGDBAPI_DIR)
if (UNIX)
  set(FileGDBAPI_HINTS_DEBUG ${FileGDBAPI_DIR}/out/x64.debug $ENV{FileGDBAPI_DIR}/out/x64.debug)
  set(FileGDBAPI_HINTS_RELEASE ${FileGDBAPI_DIR}/out/x64.release $ENV{FileGDBAPI_DIR}/out/x64.release)
  set(FileGDBAPI_PATHS /opt/libFileGDBAPI-$ENV{FileGDBAPI_VER})
elseif(WIN32)
  set(FileGDBAPI_HINTS_DEBUG ${FileGDBAPI_DIR}/build/Debug $ENV{FileGDBAPI_DIR}/build/Debug)
  set(FileGDBAPI_HINTS_RELEASE ${FileGDBAPI_DIR}/out/build/Release $ENV{FileGDBAPI_DIR}/build/Release)
else()
endif()

find_library(FileGDBAPI_LIBRARY_RELEASE
  NAMES ${FileGDBAPI_NAMES_RELEASE}
  PATH_SUFFIXES lib
  HINTS ${FileGDBAPI_HINTS} ${FileGDBAPI_HINTS_RELEASE}
  PATHS ${FileGDBAPI_PATHS}
  DOC "Google FileGDBAPI JavaScript Engine Library (Release)")

find_library(FileGDBAPI_LIBRARY_DEBUG
  NAMES ${FileGDBAPI_NAMES_DEBUG}
  PATH_SUFFIXES lib
  HINTS ${FileGDBAPI_HINTS} ${FileGDBAPI_HINTS_DEBUG}
  PATHS ${FileGDBAPI_PATHS}
  DOC "Google FileGDBAPI JavaScript Engine Library (Debug)")

if(CMAKE_BUILD_TYPE EQUAL "Release")
  get_filename_component(FileGDBAPI_LIBRARY_DIR ${FileGDBAPI_LIBRARY_RELEASE} PATH)
  set(FileGDBAPI_LIBRARY ${FileGDBAPI_LIBRARY_RELEASE})
else()
  get_filename_component(FileGDBAPI_LIBRARY_DIR ${FileGDBAPI_LIBRARY_DEBUG} PATH)
  set(FileGDBAPI_LIBRARY_DEBUG ${FileGDBAPI_LIBRARY_RELEASE})
endif()

set(FileGDBAPI_LIBRARY
  optimized ${FileGDBAPI_LIBRARY_RELEASE} debug ${FileGDBAPI_LIBRARY_DEBUG})

include(FindPackageHandleStandardArgs)

# handle the QUIETLY and REQUIRED arguments and set FileGDBAPI_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(FileGDBAPI DEFAULT_MSG
                                  FileGDBAPI_LIBRARY
				  FileGDBAPI_INCLUDE_DIR)

# Detect FileGDBAPI version
set(FileGDBAPI_VERSION_MAJOR "1")
set(FileGDBAPI_VERSION_MINOR "5")
set(FileGDBAPI_VERSION_PATCH "2")
set(FileGDBAPI_VERSION_TWEAK "")
set(FileGDBAPI_VERSION "${FileGDBAPI_VERSION_MAJOR}.${FileGDBAPI_VERSION_MINOR}.${FileGDBAPI_VERSION_PATCH}.${FileGDBAPI_VERSION_TWEAK}")

set(FileGDBAPI_VERSION_HEX 0x0${FileGDBAPI_VERSION_MAJOR}${FileGDBAPI_VERSION_MINOR}${FileGDBAPI_VERSION_PATCH}${FileGDBAPI_VERSION_TWEAK})
string(LENGTH "${FileGDBAPI_VERSION_HEX}" FileGDBAPI_VERSION_HEX_LENGTH)

while(FileGDBAPI_VERSION_HEX_LENGTH LESS 8)
  set(FileGDBAPI_VERSION_HEX "${FileGDBAPI_VERSION_HEX}0")
  string(LENGTH "${FileGDBAPI_VERSION_HEX}" FileGDBAPI_VERSION_HEX_LENGTH)
endwhile()

mark_as_advanced(FileGDBAPI_INCLUDE_DIR FileGDBAPI_LIBRARY FileGDBAPI_LIBRARY_DIR)

if(FileGDBAPI_CMAKE_DEBUG)
  message(STATUS "FileGDBAPI_INCLUDE_DIR: ${FileGDBAPI_INCLUDE_DIR}")
  message(STATUS "FileGDBAPI_LIBRARY: ${FileGDBAPI_LIBRARY}")
  message(STATUS "FileGDBAPI_LIBRARY_DEPENDS: ${FileGDBAPI_LIBRARY_DEPENDS}")
  message(STATUS "FileGDBAPI_VERSION: ${FileGDBAPI_VERSION}")
  message(STATUS "FileGDBAPI_VERSION_HEX: ${FileGDBAPI_VERSION_HEX}")
endif()
