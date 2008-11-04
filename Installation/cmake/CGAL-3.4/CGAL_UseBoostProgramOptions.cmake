if ( NOT Boost_PROGRAM_OPTIONS_FOUND )
  
  if ( NOT BUILD_SHARED_LIBS )
    set(Boost_USE_STATIC_LIBS ON)
  endif()
  
  set(Boost_FIND_VERSION 1.33.1 )
  set(Boost_FIND_VERSION_MAJOR 1 )
  set(Boost_FIND_VERSION_MINOR 33 )
  set(Boost_FIND_VERSION_PATCH 1 )
  
  find_package( Boost QUIET COMPONENTS program_options )
  
  if ( Boost_PROGRAM_OPTIONS_FOUND )
  
    if( CGAL_AUTO_LINK_ENABLED )
      message( STATUS "Boost.ProgramOptions library: found" )
    else()
      message( STATUS "Boost.ProgramOptions library: ${Boost_PROGRAM_OPTIONS_LIBRARY}" )
    endif()
        
    add_definitions( "-DCGAL_USE_BOOST_PROGRAM_OPTIONS" )
    
    if ( NOT CGAL_AUTO_LINK_ENABLED )
      link_libraries( ${Boost_PROGRAM_OPTIONS_LIBRARY} )
    endif()
  
  else()
  
    message( STATUS "Could not find Boost.ProgramOptions library" )
    
  endif()
  
endif()

