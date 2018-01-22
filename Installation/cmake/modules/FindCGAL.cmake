#
# The following module is based on FindVTK.cmake
#

# - Find a CGAL installation or binary tree.
# The following variables are set if CGAL is found.  If CGAL is not
# found, CGAL_FOUND is set to false.
#
#  CGAL_FOUND         - Set to true when CGAL is found.
#  CGAL_USE_FILE      - CMake file to use CGAL.
#

# Construct consitent error messages for use below.
set(CGAL_DIR_DESCRIPTION "directory containing CGALConfig.cmake. This is either the binary directory where CGAL was configured or PREFIX/lib/CGAL for an installation.")
set(CGAL_DIR_MESSAGE     "CGAL not found.  Set the CGAL_DIR cmake variable or environment variable to the ${CGAL_DIR_DESCRIPTION}")
 
if ( NOT CGAL_DIR )
  
  # Get the system search path as a list.
  if(UNIX)
    string(REGEX MATCHALL "[^:]+" CGAL_DIR_SEARCH1 "$ENV{PATH}")
  else()
    string(REGEX REPLACE "\\\\" "/" CGAL_DIR_SEARCH1 "$ENV{PATH}")
  endif()
  
  string(REGEX REPLACE "/;" ";" CGAL_DIR_SEARCH2 "${CGAL_DIR_SEARCH1}")

  # Construct a set of paths relative to the system search path.
  set(CGAL_DIR_SEARCH "")
  
  foreach(dir ${CGAL_DIR_SEARCH2})
  
    set(CGAL_DIR_SEARCH ${CGAL_DIR_SEARCH} ${dir}/../lib/CGAL )
      
  endforeach()


  #
  # Look for an installation or build tree.
  #
  find_path(CGAL_DIR CGALConfig.cmake

    # Look for an environment variable CGAL_DIR.
    $ENV{CGAL_DIR}

    # Look in places relative to the system executable search path.
    ${CGAL_DIR_SEARCH}

    # Look in standard UNIX install locations.
    /usr/local/lib/CGAL
    /usr/lib/CGAL

    # Read from the CMakeSetup registry entries.  It is likely that
    # CGAL will have been recently built.
    [HKEY_CURRENT_USER\\Software\\Kitware\\CMakeSetup\\Settings\\StartPath;WhereBuild1]
    [HKEY_CURRENT_USER\\Software\\Kitware\\CMakeSetup\\Settings\\StartPath;WhereBuild2]
    [HKEY_CURRENT_USER\\Software\\Kitware\\CMakeSetup\\Settings\\StartPath;WhereBuild3]
    [HKEY_CURRENT_USER\\Software\\Kitware\\CMakeSetup\\Settings\\StartPath;WhereBuild4]
    [HKEY_CURRENT_USER\\Software\\Kitware\\CMakeSetup\\Settings\\StartPath;WhereBuild5]
    [HKEY_CURRENT_USER\\Software\\Kitware\\CMakeSetup\\Settings\\StartPath;WhereBuild6]
    [HKEY_CURRENT_USER\\Software\\Kitware\\CMakeSetup\\Settings\\StartPath;WhereBuild7]
    [HKEY_CURRENT_USER\\Software\\Kitware\\CMakeSetup\\Settings\\StartPath;WhereBuild8]
    [HKEY_CURRENT_USER\\Software\\Kitware\\CMakeSetup\\Settings\\StartPath;WhereBuild9]
    [HKEY_CURRENT_USER\\Software\\Kitware\\CMakeSetup\\Settings\\StartPath;WhereBuild10]

    # Help the user find it if we cannot.
    DOC "The ${CGAL_DIR_DESCRIPTION}"
  )
  
endif()

if ( CGAL_DIR )
  
  if ( EXISTS "${CGAL_DIR}/CGALConfig.cmake" )
    include( "${CGAL_DIR}/CGALConfig.cmake" )
    set( CGAL_FOUND TRUE )
  endif()

endif()

if( NOT CGAL_FOUND)
  if(CGAL_FIND_REQUIRED)
    MESSAGE(FATAL_ERROR ${CGAL_DIR_MESSAGE})
  else()
    if(NOT CGAL_FIND_QUIETLY)
      MESSAGE(STATUS ${CGAL_DIR_MESSAGE})
    endif()
  endif()
endif()
