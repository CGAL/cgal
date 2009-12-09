if( NOT CGAL_RS_SETUP )

  option( WITH_RS "Search for RS" ON )
  #option( WITH_RS3 "Search for RS3" OFF )

  if( WITH_RS )

    find_package( MPFI )
    find_package( RS )

    if( RS_FOUND AND MPFI_FOUND )

      include(CGAL_Macros)

      if(CMAKE_OSX_ARCHITECTURES STREQUAL "ppc")
        message( STATUS "TLS is not supported on this architecture" )
        add_to_cached_list(CGAL_3RD_PARTY_DEFINITIONS "-DCGAL_RS_NO_TLS" )
      endif(CMAKE_OSX_ARCHITECTURES STREQUAL "ppc")

      set ( CGAL_USE_RS 1 )

      # add rs3 parameters, if necessary (rs3 must be always before rsexport)
      if(WITH_RS3 AND RS3_FOUND)
        if(CMAKE_OSX_ARCHITECTURES STREQUAL "ppc")
          message( STATUS "RS3 is not supported on this architecture" )
        else(CMAKE_OSX_ARCHITECTURES STREQUAL "ppc")
          add_to_cached_list(CGAL_3RD_PARTY_INCLUDE_DIRS ${RS3_INCLUDE_DIR})
          add_to_cached_list(CGAL_3RD_PARTY_LIBRARIES ${RS3_LIBRARIES})
          add_to_cached_list(CGAL_3RD_PARTY_DEFINITIONS "-DCGAL_USE_RS3")
        endif(CMAKE_OSX_ARCHITECTURES STREQUAL "ppc")
      endif()

      # add rsexport parameters
      add_to_cached_list(CGAL_3RD_PARTY_INCLUDE_DIRS ${RS_INCLUDE_DIR})
      add_to_cached_list(CGAL_3RD_PARTY_LIBRARIES ${RS_LIBRARIES})
      add_to_cached_list(CGAL_3RD_PARTY_DEFINITIONS ${RS_DEFINITIONS})

      uniquely_add_flags(CMAKE_CXX_FLAGS ${RS_CXX_FLAGS})

      if(BUILD_SHARED_LIBS)
        uniquely_add_flags(CMAKE_SHARED_LINKER_FLAGS ${RS_LINKER_FLAGS})
      else(BUILD_SHARED_LIBS)
        uniquely_add_flags(CMAKE_MODULE_LINKER_FLAGS ${RS_LINKER_FLAGS})
      endif(BUILD_SHARED_LIBS)

    endif( RS_FOUND AND MPFI_FOUND )

  endif( WITH_RS )

  set( CGAL_RS_SETUP TRUE )

endif( NOT CGAL_RS_SETUP )
