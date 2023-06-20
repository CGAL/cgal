include(${CMAKE_CURRENT_LIST_DIR}/CGAL_Macros.cmake)

option(CGAL_DEV_MODE
  "Activate the CGAL developers mode. See https://github.com/CGAL/cgal/wiki/CGAL_DEV_MODE"
  $ENV{CGAL_DEV_MODE})

if(RUNNING_CGAL_AUTO_TEST OR CGAL_TEST_SUITE)
# Just to avoid a warning from CMake if that variable is set on the command line...
endif()

# Common settings for CGAL cmake scripts
if( NOT CGAL_COMMON_FILE_INCLUDED )
  set(CGAL_COMMON_FILE_INCLUDED 1 )

  # CMAKE_ROOT must be properly configured, but is not by the CMake windows installer, so check here
  if (NOT CMAKE_ROOT)
    message( FATAL_ERROR "CMAKE_ROOT environment variable not set. It should point to the directory where CMake is installed.")
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
