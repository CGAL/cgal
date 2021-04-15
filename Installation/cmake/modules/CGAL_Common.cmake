include(${CMAKE_CURRENT_LIST_DIR}/CGAL_Macros.cmake)

option(CGAL_DEV_MODE
  "Activate the CGAL developers mode. See https://github.com/CGAL/cgal/wiki/CGAL_DEV_MODE"
  $ENV{CGAL_DEV_MODE})

if(RUNNING_CGAL_AUTO_TEST)
# Just to avoid a warning from CMake if that variable is set on the command line...
endif()

# Common settings for CGAL cmake scripts
if( NOT CGAL_COMMON_FILE_INCLUDED )
  set(CGAL_COMMON_FILE_INCLUDED 1 )

  # CMAKE_ROOT must be properly configured, but is not by the CMake windows installer, so check here
  if (NOT CMAKE_ROOT)
    message( FATAL_ERROR "CMAKE_ROOT environment variable not set. It should point to the directory where CMake is installed.")
  endif()

  # CMAKE_VERSION was introduced in 2.6.3 so we use it to detect the fact
  if ( CMAKE_VERSION )
    set( CMAKE_2_6_3_OR_ABOVE TRUE )
  else()
    set( CMAKE_2_6_3_OR_ABOVE FALSE )
  endif()

  if ( CGAL_BUILDING_LIBS )
    option(BUILD_SHARED_LIBS "Build shared libraries" ON)
    set(CGAL_BUILD_SHARED_LIBS ${BUILD_SHARED_LIBS})

    if ( BUILD_SHARED_LIBS )
      message( STATUS "Building shared libraries" )
    else()
      message( STATUS "Building static libraries" )
    endif()
  endif()

  if ( WIN32 )
    find_program(CMAKE_UNAME uname /bin /usr/bin /usr/local/bin )
    if(CMAKE_UNAME)
      exec_program(uname ARGS -s OUTPUT_VARIABLE CMAKE_SYSTEM_NAME2)
      if ( CMAKE_SYSTEM_NAME2 MATCHES "CYGWIN" )
        message( STATUS "This is the Windows CMake running within the cygwin platform." )
        set( CGAL_WIN32_CMAKE_ON_CYGWIN TRUE CACHE INTERNAL "This is the cygwin platform." )
      endif()
    endif()
    hide_variable(CMAKE_UNAME)

    # Optionally setup the Visual Leak Detector
    include(${CMAKE_CURRENT_LIST_DIR}/CGAL_SetupVLD.cmake)
    CGAL_SetupVLD()
    if(VLD_FOUND)
      message(STATUS "Visual Leak Detector (VLD) is enabled.")
    else()
      message(STATUS "Visual Leak Detector (VLD) is not found.")
    endif()
  endif()

  # set minimal version of some optional libraries:
  set( Eigen3_FIND_VERSION "3.1.0")
  # set use-file for Eigen3 (needed to have default solvers)
  set(EIGEN3_USE_FILE "UseEigen3")

endif()
