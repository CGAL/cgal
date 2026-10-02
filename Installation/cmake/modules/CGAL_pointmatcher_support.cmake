if(libpointmatcher_FOUND AND NOT TARGET CGAL::pointmatcher_support)
  if (libpointmatcher_VERSION VERSION_GREATER_EQUAL "1.4.4")
    find_package(Boost COMPONENTS thread program_options date_time chrono)
   else()
    find_package(Boost COMPONENTS thread filesystem system program_options date_time chrono)
  endif()

  if(Boost_chrono_FOUND
       AND Boost_thread_FOUND
       AND Boost_program_options_FOUND
       AND Boost_date_time_FOUND
       AND (libpointmatcher_VERSION VERSION_GREATER_EQUAL "1.4.4" OR (Boost_filesystem_FOUND AND Boost_system_FOUND)))
    add_library(CGAL::pointmatcher_support INTERFACE IMPORTED)
    target_compile_options(CGAL::pointmatcher_support INTERFACE "-D_USE_MATH_DEFINES")
    target_compile_definitions(CGAL::pointmatcher_support INTERFACE "CGAL_LINKED_WITH_POINTMATCHER")
    target_include_directories(CGAL::pointmatcher_support INTERFACE "${libpointmatcher_INCLUDE_DIRS}")
    target_link_libraries(CGAL::pointmatcher_support INTERFACE ${libpointmatcher_LIBRARIES} libnabo::nabo)
  else()
    if (libpointmatcher_VERSION VERSION_GREATER_EQUAL "1.4.4")
      message(STATUS "NOTICE: the libpointmatcher library requires the following boost components: thread program_options date_time chrono.")
    else()
      message(STATUS "NOTICE: the libpointmatcher library requires the following boost components: thread filesystem system program_options date_time chrono.")
    endif()
  endif()
endif()
