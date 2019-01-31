pushd AABB_tree_Demo;                   zip -r ../AABB_demo.zip          *    ;      popd
pushd Alpha_shapes_3_Demo;              zip -r ../alpha_shape_3.zip      *    ;      popd
pushd Bounding_volumes_Demo;            zip -r ../bounding_volumes_2.zip *    ;      popd
pushd Circular_kernel_2_Demo;           zip -r ../circular_kernel.zip    *    ;      popd
pushd Mesh_3_Demo;                      zip -r ../mesh_3.zip             *    ;      popd
pushd Periodic_3_triangulation_3_Demo;  zip -r ../periodic_3_triangulation_3.zip *;  popd
pushd Periodic_Lloyd_3_Demo;            zip -r ../periodic_3_triangulation_3.zip *;  popd
pushd Polygon_Demo;                     zip -r ../polygon.zip *               ;      popd
pushd Polyhedron_Demo;                  zip -r ../polyhedron_3.zip *          ;      popd
pushd Segment_Delaunay_graph_2_Demo;    zip -r ../segment_voronoi_diagram_2.zip *;   popd
pushd Surface_mesher_Demo;              zip -r ../surface_mesher.zip *;              popd
pushd Triangulation_2_Demo;
  zip -r ../regular_triangulation_2.zip Regular_triangulation_2.exe
  zip -r ../constrained_delaunay_triangulation_2.zip Constrained_Delaunay_triangulation_2.exe
  zip -r ../delaunay_triangulation_2.zip Delaunay_triangulation_2.exe
popd
