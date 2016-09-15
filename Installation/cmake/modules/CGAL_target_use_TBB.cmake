if (CGAL_target_use_TBB_included)
  return()
endif()
set(CGAL_target_use_TBB_included TRUE)

function(CGAL_target_use_TBB target)
  set(keyword PUBLIC)

  target_include_directories ( ${target} SYSTEM ${keyword} ${TBB_INCLUDE_DIRS} )
  target_link_libraries( ${target} ${keyword} ${TBB_LIBRARIES} )
  target_compile_options( ${target} ${keyword} -DNOMINMAX -DCGAL_LINKED_WITH_TBB )
endfunction()
