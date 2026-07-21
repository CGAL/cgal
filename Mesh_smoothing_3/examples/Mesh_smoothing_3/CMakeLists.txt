

create_single_source_cgal_program("cgal_offset_computation.cpp")
create_single_source_cgal_program("cgal_c3t3_test.cpp")
target_link_libraries(cgal_offset_computation PRIVATE Mesh_smoothing_3)
target_link_libraries(cgal_c3t3_test PRIVATE Mesh_smoothing_3)



