if(TBB_FOUND AND NOT TARGET CGAL::TBB_support)
  if(NOT TARGET Threads::Threads)
    find_package(Threads REQUIRED)
  endif()
  add_library(CGAL::TBB_support INTERFACE IMPORTED)
  set_target_properties(CGAL::TBB_support PROPERTIES
    INTERFACE_COMPILE_DEFINITIONS "CGAL_LINKED_WITH_TBB;NOMINMAX"
    INTERFACE_INCLUDE_DIRECTORIES "${TBB_INCLUDE_DIRS}"
    INTERFACE_LINK_LIBRARIES "TBB::tbb;TBB::tbbmalloc;Threads::Threads")
    if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.25)
      set_target_properties(CGAL::TBB_support PROPERTIES SYSTEM TRUE)
    endif()
endif()
