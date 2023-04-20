The CGAL Open Source Project is pleased to announce the release 5.5 Beta 1 of CGAL, the Computational Geometry Algorithms Library.

CGAL version 5.5 Beta 1 is a public testing release. It should provide a solid ground to report bugs that need to be tackled before the release of the final version of CGAL 5.5 in July 2022.

Besides fixes and general enhancement to existing packages, the following has changed since CGAL 5.4:

### [3D Alpha Wrapping (new package)](https://doc.cgal.org/5.5/Manual/packages.html#PkgAlphaWrap3)

-   This component takes a 3D triangle mesh, soup, or point set as input, and generates a valid (watertight, intersection-free, and combinatorially 2-manifold) surface triangle mesh that contains the input.
    The algorithm proceeds by shrink-wrapping and refining a 3D Delaunay triangulation, starting from a loose bounding box of the input.
    Two user-defined parameters, alpha and offset, offer control over the maximum size of cavities where the shrink-wrapping process can enter, and the tightness of the final surface mesh to the input, respectively. Once combined, these parameters provide a means to trade fidelity
    to the input for complexity of the output.

    See also the [announcement page](https://www.cgal.org/2022/05/18/alpha_wrap/).

### [3D Convex Hulls](https://doc.cgal.org/5.5/Manual/packages.html#PkgConvexHull3)

-   Added an [overload of the function `CGAL::convex_hull_3()`](https://doc.cgal.org/5.5/Convex_hull_3/group__PkgConvexHull3Functions.html#ga52fca4745c2ef0351063fbe66b035fd1), which writes the result in an indexed triangle set.

### [2D Polygons](https://doc.cgal.org/5.5/Manual/packages.html#PkgPolygon2)

-   Add vertex, edge, and hole ranges.
-   The concept [`GeneralPolygonWithHoles_2`](https://doc.cgal.org/5.5/Polygon/classGeneralPolygonWithHoles__2.html) now requires the nested type `Polygon_2` instead of `General_polygon_2`.

### [2D Regularized Boolean Set-Operations](https://doc.cgal.org/5.5/Manual/packages.html#PkgBooleanSetOperations2)
-   The concept [`GeneralPolygonSetTraits_2`](https://doc.cgal.org/5.5/Boolean_set_operations_2/classGeneralPolygonSetTraits__2.html) now requires the nested type `Construct_polygon_with_holes_2` instead of `Construct_general_polygon_with_holes_2`.

### [Combinatorial Maps](https://doc.cgal.org/5.5/Manual/packages.html#PkgCombinatorialMaps)

- Removed old code deprecated in CGAL 4.9 and 4.10 (global functions, and information associated with darts).

### [2D Arrangements](https://doc.cgal.org/5.5/Manual/packages.html#PkgArrangementOnSurface2)
-   Fixed the `intersect_2`, `compare_y_at_x_right`, and `compare_y_at_x_left` function objects of the traits class template [`Arr_geodesic_arc_on_sphere_traits_2`](https://doc.cgal.org/5.5/Arrangement_on_surface_2/classCGAL_1_1Arr__geodesic__arc__on__sphere__traits__2.html) that handles geodesic arcs on sphere and applied a small syntactical fix to the tracing traits.

### [Tetrahedral Mesh Generation](https://doc.cgal.org/5.5/Manual/packages.html#PkgMesh3)

-   Added the function [`remove_isolated_vertices()`](https://doc.cgal.org/5.5/Mesh_3/classCGAL_1_1Mesh__complex__3__in__triangulation__3.html#ace57c4e777da457c6e33b4f6e89949ce) as a post-processing step for the tetrahedral mesh generation.

### [Polygon Mesh Processing](https://doc.cgal.org/5.5/Manual/packages.html#PkgPolygonMeshProcessing)
-   Added the function [`CGAL::Polygon_mesh_processing::orient_triangle_soup_with_reference_triangle_soup()`](https://doc.cgal.org/5.5/Polygon_mesh_processing/group__PMP__orientation__grp.html#ga855b1c55c201b91ab04eebd2811a87fd), which enables re-orienting the faces of a triangle soup based on the orientation of the nearest face in a reference triangle soup.
-   Added the function [`CGAL::Polygon_mesh_processing::compatible_orientations()`](https://doc.cgal.org/5.5/Polygon_mesh_processing/group__PMP__orientation__grp.html#ga9ac9b9434084b64f3304df636c3178a3), which enables to retrieve the (in)compatibility of orientations of faces from different connected components.
-   Added the function [`CGAL::Polygon_mesh_processing::tangential_relaxation()`](https://doc.cgal.org/5.5/Polygon_mesh_processing/group__PMP__meshing__grp.html#ga136c659162e5360354db5879db7431b4), which applies an area-based tangential mesh smoothing to the vertices of a surface triangle mesh.
-   Added the named parameter `visitor` to the function [`triangulate_hole()`](https://doc.cgal.org/5.5/Polygon_mesh_processing/group__PMP__hole__filling__grp.html#gad2d3c43bce0ef90a16530478196d7f42), which enables to track progress with callbacks.
-   Added more functions in the [visitor of the corefinement based methods](https://doc.cgal.org/5.5/Polygon_mesh_processing/classPMPCorefinementVisitor.html) to track progress.

### [Surface Mesh Simplification](https://doc.cgal.org/5.5/Manual/packages.html#PkgSurfaceMeshSimplification)
-   Introduced four variations of the Garland-Heckbert simplification algorithm based on the probabilistic approach of Trettner and Kobbelt (Fast and Robust QEF Minimization using Probabilistic Quadrics): [`GarlandHeckbert_plane_policies`](https://doc.cgal.org/5.5/Surface_mesh_simplification/classCGAL_1_1Surface__mesh__simplification_1_1GarlandHeckbert__plane__policies.html), [`GarlandHeckbert_probabilistic_plane_policies`](https://doc.cgal.org/5.5/Surface_mesh_simplification/classCGAL_1_1Surface__mesh__simplification_1_1GarlandHeckbert__probabilistic__plane__policies.html), [`GarlandHeckbert_triangle_policies`](https://doc.cgal.org/5.5/Surface_mesh_simplification/classCGAL_1_1Surface__mesh__simplification_1_1GarlandHeckbert__triangle__policies.html), and [`GarlandHeckbert_probabilistic_triangle_policies`](https://doc.cgal.org/5.5/Surface_mesh_simplification/classCGAL_1_1Surface__mesh__simplification_1_1GarlandHeckbert__probabilistic__triangle__policies.html).
-   The class `GarlandHeckbert_policies` has been deprecated, `GarlandHeckbert_plane_policies` replaces it.

### [Point Set Processing](https://doc.cgal.org/5.5/Manual/packages.html#PkgPointSetProcessing3)

-   A new optional named parameter, `min_points_per_cell` has been added to [`grid_simplify_point_set()`](https://doc.cgal.org/5.5/Point_set_processing_3/group__PkgPointSetProcessing3Algorithms.html#ga7757ef9b3900e42fde26f5a0ac56e20f). By adding a minimal number of points in a cell such that a point is retained, one can also filter out low density areas and outliers: in the case of densely sampled point clouds, this yields better results than using grid simplification and then outlier removal, while being very vast. The default value is `1` to keep the previous behavior as default.

### [dD Spatial Searching](https://doc.cgal.org/5.5/Manual/packages.html#PkgSpatialSearchingD)

-   Added the member function [`write_graphviz()`](https://doc.cgal.org/5.5/Spatial_searching/classCGAL_1_1Kd__tree.html#ac2851b5cafb8d5cce0dc5fb107c8f13f) to the class `Kd_tree` that writes the tree in a stream in the [Graphviz](https://graphviz.org/) format.

### [CGAL and the Boost Graph Library (BGL)](https://doc.cgal.org/5.5/Manual/packages.html#PkgBGL)

-   Added the function [`invert_selection()`](https://doc.cgal.org/5.5/BGL/structCGAL_1_1Face__filtered__graph.html#aa428541ebbdd35f9a6e9a3ffd60178df) in the class [`Face_filtered_graph`](https://doc.cgal.org/5.5/BGL/structCGAL_1_1Face__filtered__graph.html), which toggles the selected status of a graph: selected faces are deselected, and unselected faces are selected.


