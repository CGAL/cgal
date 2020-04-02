if (CGAL_target_use_GLPK_included)
  return()
endif()
set(CGAL_target_use_GLPK_included TRUE)

function(CGAL_target_use_GLPK target)
  target_include_directories(${target} PUBLIC ${GLPK_INCLUDE_DIR})
  target_compile_options(${target} PUBLIC -DCGAL_USE_GLPK)
  target_link_libraries(${target} PUBLIC ${GLPK_LIBRARIES})
endfunction()

