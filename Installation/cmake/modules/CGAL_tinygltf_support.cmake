if(TINYGLTF_FOUND AND NOT TARGET CGAL::tinygltf_support)
  add_library(CGAL::tinygltf_support INTERFACE IMPORTED)
  #definition only needed for v2 series, off by default in v3
  target_compile_definitions(CGAL::tinygltf_support INTERFACE TINYGLTF_NO_STB_IMAGE)
  target_compile_definitions(CGAL::tinygltf_support INTERFACE TINYGLTF_NO_STB_IMAGE_WRITE)
  set_target_properties(CGAL::tinygltf_support PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES "${TINYGLTF_INCLUDE_DIR}")
endif()