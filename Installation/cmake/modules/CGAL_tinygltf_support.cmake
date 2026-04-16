if(TINYGLTF_FOUND AND NOT TARGET CGAL::tinygltf_support)
  add_library(CGAL::tinygltf_support INTERFACE IMPORTED)
  set_target_properties(CGAL::tinygltf_support PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES "${TINYGLTF_INCLUDE_DIR}")
endif()