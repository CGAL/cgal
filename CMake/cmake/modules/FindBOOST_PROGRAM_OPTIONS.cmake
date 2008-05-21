# Try to find the BOOST_PROGRAM_OPTIONS library
# BOOST_PROGRAM_OPTIONS_FOUND	    - system has BOOST_PROGRAM_OPTIONS lib
# BOOST_PROGRAM_OPTIONS_LIBRARIES   - libraries needed to use BOOST_PROGRAM_OPTIONS

# TODO: support MacOSX

# BOOST_PROGRAM_OPTIONS needs Boost
include(FindPackageHandleStandardArgs)

if(Boost_FOUND AND Boost_LIBRARY_DIRS)
  if (BOOST_PROGRAM_OPTIONS_LIBRARIES)
      # Already in cache, be silent
      set(BOOST_PROGRAM_OPTIONS_FIND_QUIETLY TRUE)
  endif()

  if ( AUTO_LINK_ENABLED )
    file ( GLOB BOOST_PROGRAM_OPTIONS_LIBRARIES "${Boost_LIBRARY_DIRS}/libboost_program_options*" )
  else()
    find_library(BOOST_PROGRAM_OPTIONS_LIBRARIES NAMES boost_program_options boost_program_options-mt PATHS ${Boost_LIBRARY_DIRS})
  endif()                                                
  
  hide_variable(BOOST_PROGRAM_OPTIONS_LIBRARIES)
  
  find_package_handle_standard_args(BOOST_PROGRAM_OPTIONS "boost_program_options not found." BOOST_PROGRAM_OPTIONS_LIBRARIES )
  
endif()

