if (CGAL_target_use_Boost_Serialization_included)
  return()
endif()
set(CGAL_target_use_Boost_Serialization_included TRUE)

function(CGAL_target_use_Boost_Serialization target)

  if(TARGET Boost::serialization)
    target_link_libraries(${target} PUBLIC Boost::serialization)
  else()
    target_link_libraries(${target} PUBLIC ${Boost_SERIALIZATION_LIBRARY})
  endif()

  target_compile_options( ${target} PUBLIC -DCGAL_LINKED_WITH_BOOST_SERIALIZATION)

endfunction()
