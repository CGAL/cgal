if(Boost_SERIALIZATION_FOUND AND NOT TARGET CGAL::Boost_serialization_support)
  if(TARGET Boost::serialization)
    set(Boost_LIB Boost::serialization)
  else()
    set(Boost_LIB  ${Boost_SERIALIZATION_LIBRARY})
  endif()

  add_library(CGAL::Boost_serialization_support INTERFACE IMPORTED)
  set_target_properties(CGAL::Boost_serialization_support PROPERTIES
    INTERFACE_COMPILE_DEFINITIONS "CGAL_LINKED_WITH_BOOST_SERIALIZATION"
    INTERFACE_LINK_LIBRARIES ${Boost_LIB})
endif()
