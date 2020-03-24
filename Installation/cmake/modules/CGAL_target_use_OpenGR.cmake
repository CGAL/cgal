if (CGAL_target_use_OpenGR_included)
  return()
endif()
set(CGAL_target_use_OpenGR_included TRUE)

function(CGAL_target_use_OpenGR target)
  target_include_directories(${target} PUBLIC ${OpenGR_INCLUDE_DIR})
  target_compile_options( ${target} PUBLIC -DCGAL_LINKED_WITH_OPENGR)
endfunction()
