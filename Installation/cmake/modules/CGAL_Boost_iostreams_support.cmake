cmake_minimum_required(VERSION 3.11...3.23)
if(Boost_IOSTREAMS_FOUND AND NOT TARGET CGAL::Boost_iostreams_support)

  if( WIN32 )
    # to avoid a warning with old cmake
    set(_Boost_BZIP2_HEADERS             "boost/iostreams/filter/bzip2.hpp")
    set(_Boost_ZLIB_HEADERS              "boost/iostreams/filter/zlib.hpp")
    find_package( Boost OPTIONAL_COMPONENTS bzip2 zlib)
    if (Boost_ZLIB_FOUND AND Boost_BZIP2_FOUND)
      set(ZLIB_LIBS ${Boost_ZLIB_LIBRARY} ${Boost_BZIP2_LIBRARY})
    else()
      message(STATUS "NOTICE: This project requires Boost ZLIB and Boost BZIP2, and will not be compiled.")
      return()
    endif()

  else()

    find_package(ZLIB QUIET)
    if(ZLIB_FOUND)
      set(ZLIB_LIBS ZLIB::ZLIB)
    else()
      message(STATUS "NOTICE: This project requires ZLIB, and will not be compiled.")
      return()
    endif()

  endif()

  if(TARGET Boost::iostreams)
    set(Boost_LIB Boost::iostreams)
  else()
    set(Boost_LIB  ${Boost_IOSTREAMS_LIBRARY})
  endif()

  add_library(CGAL::Boost_iostreams_support INTERFACE IMPORTED)

  set_target_properties(CGAL::Boost_iostreams_support PROPERTIES
    INTERFACE_COMPILE_DEFINITIONS "CGAL_LINKED_WITH_BOOST_IOSTREAMS")
  target_link_libraries(CGAL::Boost_iostreams_support INTERFACE ${Boost_LIB} ${ZLIB_LIBS})
endif()
