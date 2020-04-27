if (CGAL_target_use_SCIP_included)
  return()
endif()
set(CGAL_target_use_SCIP_included TRUE)

function(CGAL_target_use_SCIP target)
  target_include_directories(${target} PUBLIC ${SCIP_INCLUDE_DIRS})
  target_compile_options(${target} PUBLIC -DCGAL_USE_SCIP)
  target_link_libraries(${target} PUBLIC ${SCIP_LIBRARIES})
endfunction()

