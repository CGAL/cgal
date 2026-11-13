if(OSQP_FOUND AND NOT TARGET CGAL::OSQP_support)
  add_library(CGAL::OSQP_support INTERFACE IMPORTED)
  set_target_properties(CGAL::OSQP_support PROPERTIES
    INTERFACE_COMPILE_DEFINITIONS "CGAL_USE_OSQP")
  if(TARGET osqp::osqp)
    target_link_libraries(CGAL::OSQP_support INTERFACE osqp::osqp)
  elseif(TARGET osqp::osqpstatic)
    target_link_libraries(CGAL::OSQP_support INTERFACE osqp::osqpstatic)
  else()
    set_target_properties(CGAL::OSQP_support PROPERTIES
      INTERFACE_INCLUDE_DIRECTORIES "${OSQP_INCLUDE_DIR}"
      INTERFACE_LINK_LIBRARIES "${OSQP_LIBRARIES}")
  endif()
endif()
