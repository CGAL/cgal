# Try to find the BOOST_PROGRAM_OPTIONS library
# BOOST_PROGRAM_OPTIONS_FOUND	    - system has BOOST_PROGRAM_OPTIONS lib
# BOOST_PROGRAM_OPTIONS_LIBRARIES   - libraries needed to use BOOST_PROGRAM_OPTIONS

# TODO: support MacOSX

# BOOST_PROGRAM_OPTIONS needs Boost
INCLUDE(FindPackageHandleStandardArgs)

find_package(Boost2 QUIET)
if(Boost_FOUND AND Boost_LIBRARY_DIRS)
  if (BOOST_PROGRAM_OPTIONS_LIBRARIES)
      # Already in cache, be silent
      set(BOOST_PROGRAM_OPTIONS_FIND_QUIETLY TRUE)
  endif()

  if ( AUTO_LINK_ENABLED )
    set(BOOST_PROGRAM_OPTIONS_LIBRARIES "")
    FIND_PACKAGE_HANDLE_STANDARD_ARGS(BOOST_PROGRAM_OPTIONS "boost_program_options not found." Boost_LIBRARY_DIRS )
  else()
    find_library(BOOST_PROGRAM_OPTIONS_LIBRARIES NAMES libboost_program_options PATHS ${Boost_LIBRARY_DIRS})
    FIND_PACKAGE_HANDLE_STANDARD_ARGS(BOOST_PROGRAM_OPTIONS "boost_program_options not found." BOOST_PROGRAM_OPTIONS_LIBRARIES )
  endif()                                                
endif()

