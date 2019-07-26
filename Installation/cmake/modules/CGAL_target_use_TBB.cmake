if (CGAL_target_use_TBB_included)
  return()
endif()
set(CGAL_target_use_TBB_included TRUE)

set(TBB_USE_FILE "UseTBB")

function(CGAL_target_use_TBB target)
  if(NOT TARGET Threads::Threads)
    find_package(Threads REQUIRED)
  endif()
  target_link_libraries( ${target} PUBLIC TBB::tbb TBB::tbbmalloc Threads::Threads)
  target_compile_options( ${target} PUBLIC -DNOMINMAX -DCGAL_LINKED_WITH_TBB )
endfunction()
