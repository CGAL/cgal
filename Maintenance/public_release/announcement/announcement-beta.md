The CGAL Open Source Project is pleased to announce the release 5.6 Beta 1 of CGAL, the Computational Geometry Algorithms Library.

CGAL version 5.6 Beta 1 is a public testing release. It should provide a solid ground to report bugs that need to be tackled before the release of the final version of CGAL 5.6 in July 2022.

Besides fixes and general enhancement to existing packages, the following has changed since CGAL 5.5:

### General Changes

-   **Breaking change**: Package-specific assertions, preconditions, and postconditions (such as `CGAL_triangulation_assertion`) have been removed. Corresponding CGAL-wide versions (such as `CGAL_assertion`) should be used instead.

### [Shape Detection](https://doc.cgal.org/5.6/Manual/packages.html#PkgShapeDetection) (major changes)

-   **Breaking change**: The region growing part of the package have been reworked to fix design issues introduced with the handling of `FaceGraph` models. In particular, the notion of `Item` has been introduced to reference an element in the input range of elements. Region maps now operates on `Item` and no longer on the value type of the input range.
-   **Breaking change**: The method `update()` in the concept `RegionType` now returns a `Boolean` instead of `void`, that is used inside the class `Region_growing` for detecting if the input conditions for the new region are satisfied. This change affects only user-defined types of regions.
-   **Breaking change**: The constructors of all models used together with the region growing algorithm now enable users to provide parameters through the [named parameters](https://doc.cgal.org/5.6/BGL/group__bgl__namedparameters.html) mechanism.
-   All fitting classes in the region growing framework are now using better versions of the region conditions, more precise and faster, including the correct normal orientations.
-   Added new models of the concept `RegionType` for getting linear regions in a set of 2D and 3D segments and on 2D and 3D polylines.
-   Added the class `Polyline_graph` for extracting a set of polylines from a face graph, which splits this graph into a set of user-defined regions.
-   Added new shapes to the Region Growing algorithm on a point set: circles in 2D, spheres in 3D, and cylinders in 3D.

### [2D Straight Skeleton and Polygon Offsetting](https://doc.cgal.org/5.6/Manual/packages.html#PkgStraightSkeleton2) (major changes)
-   Added weighted straight skeletons: weighted straight skeletons are a generalization of straight skeletons. Contour edges are assigned a positive weight, which can be understood as assigning a speed to the wavefront spawned from the contour edge.
-   Added straight skeleton extrusion: this CGAL package now implements the extrusion of weighted straight skeletons of polygons with holes. The output is a closed, combinatorially 2-manifold surface triangle mesh.
 See also the [news entry](https://www.cgal.org/2023/05/09/improved_straight_skeleton/).

### [2D and 3D Linear Geometry Kernel](https://doc.cgal.org/5.6/Manual/packages.html#PkgKernel23)

-   Added the functor [`CompareAngle_3`](https://doc.cgal.org/5.6/Kernel_23/classKernel_1_1CompareAngle__3.html) to the concept [`Kernel`](https://doc.cgal.org/5.6/Kernel_23/classKernel.html) to compare an angle defined by three points to the cosinus of another angle.

### [Combinatorial Maps](https://doc.cgal.org/5.6/Manual/packages.html#PkgCombinatorialMaps), [Generalized Maps](https://doc.cgal.org/5.6/Manual/packages.html#PkgGeneralizedMaps), and [Linear Cell Complex](https://doc.cgal.org/5.6/Manual/packages.html#PkgLinearCellComplex)

-   Added a version that uses indices instead of handles as dart and attribute descriptors. As the indices are integers convertible from and to `std::size_t`, they can be used as index into vectors which store properties. To use the index version, `Use_index` must be defined and be equal to `CGAL::Tag_true` in the item class.

### [Linear Cell Complex](https://doc.cgal.org/5.6/Manual/packages.html#PkgLinearCellComplex)

-   Added the class [`Linear_cell_complex_incremental_builder_3`](https://doc.cgal.org/5.6/Linear_cell_complex/classCGAL_1_1Linear__cell__complex__incremental__builder__3.html).

### [2D Arrangements](https://doc.cgal.org/5.6/Manual/packages.html#PkgArrangementOnSurface2)

-   Introduced an overload function template, namely `draw(arr)`, that renders arrangements based on the `Basic_viewer_qt` class template. As of now, only 2D arrangements on the plane induced by (i) segments, (ii) conics, and (iii) circular arcs or (linear) segments are supported.
-   Improved the traits class template that handles conics, namely [`Arr_conic_traits_2`](https://doc.cgal.org/5.6/Arrangement_on_surface_2/classCGAL_1_1Arr__conic__traits__2.html). This includes the following: 1. Fixed a couple of bugs and slightly optimized some functions. 2. Introduced functionality that approximates conics with polylines. (This is used to draw conic curves.) 3. **Breaking change**: Changed the interface to generate conic curves. In the past, curves where generated directly using the constructors of the conic and x-monotone conic constructs. Now, they are constructed via function objects provided by the traits. This eliminates the constructions of temporary kernels. The old functionality is obsolete, but still supported for a limited number of versions. It depends on a static member function of the traits. In a future version this function will no longer be static, implying that the old functionality will no longer be supported.
- Introduced functionality that approximates circular segments with polylines. (This is used to draw conic curves.)

### [Polygon Mesh Processing](https://doc.cgal.org/5.6/Manual/packages.html#PkgPolygonMeshProcessing)

-   Added functions [`CGAL::Polygon_mesh_processing::region_growing_of_planes_on_faces()`](https://doc.cgal.org/5.6/Polygon_mesh_processing/group__PkgPolygonMeshProcessingRef.html#ga50dcd2f6295f584d2e378b57290ae2af) and [`CGAL::Polygon_mesh_processing::detect_corners_of_regions()`](https://doc.cgal.org/5.6/Polygon_mesh_processing/group__PkgPolygonMeshProcessingRef.html#gac8e445730d718a2fc49604e865017d2e), which enable partitioning a mesh into planar regions using the region growing algorithm from the [Shape Detection](https://doc.cgal.org/5.6/Manual/packages.html#PkgShapeDetection) package.

-   Added the functions [`CGAL::Polygon_mesh_processing::remesh_planar_patches()`](https://doc.cgal.org/5.6/Polygon_mesh_processing/group__PMP__meshing__grp.html#ga7fca6fa2db94560ab6d32e6a77fc35b6) and [`CGAL::Polygon_mesh_processing::remesh_almost_planar_patches()`](https://doc.cgal.org/5.6/Polygon_mesh_processing/group__PMP__meshing__grp.html#ga0e6da479548199a5d82c3cf0ed36e8a0), which can be used to remesh patches of coplanar faces in a mesh.

-   Added the function [`CGAL::Polygon_mesh_processing::surface_Delaunay_remeshing()`](https://doc.cgal.org/5.6/Polygon_mesh_processing/group__PMP__meshing__grp.html#gaff62f9415d2fe96d1d3095351f156ced), which can be used to remesh a surface triangle mesh using the Delaunay refinement algorithm from the [3D Mesh Generation](https://doc.cgal.org/5.6/Manual/packages.html#PkgMesh3) package.

-   Added the function [`CGAL::Polygon_mesh_processing::remove_almost_degenerate_faces()`](https://doc.cgal.org/5.6/Polygon_mesh_processing/group__PMP__geometric__repair__grp.html#ga48008d2b66de8a68a7068f29db15dad6), which can be used to remove badly shaped triangles faces in a mesh.

-   Added the functions [`CGAL::Polygon_mesh_processing::does_triangle_soup_self_intersect()`](https://doc.cgal.org/5.6/Polygon_mesh_processing/group__PMP__intersection__grp.html#ga4909920dc48b8285e69feb845feb1e53) and [`CGAL::Polygon_mesh_processing::triangle_soup_self_intersections()`](https://doc.cgal.org/5.6/Polygon_mesh_processing/group__PMP__intersection__grp.html#ga1c5fee17bd0d92d5a2fba77ed94d4b4d) to identify and report self-intersections in a triangle soup, similarly to existing functions on triangle meshes.

-   Added the function [`CGAL::Polygon_mesh_processing::triangulate_polygons()`](https://doc.cgal.org/5.6/Polygon_mesh_processing/group__PMP__meshing__grp.html#ga8b7db6aa8c3e79526b594739ba926d82), which allows users to triangulate polygon soups.

-   Added a named parameter to [`CGAL::Polygon_mesh_processing::smooth_shape()`](https://doc.cgal.org/5.6/Polygon_mesh_processing/group__PMP__meshing__grp.html#ga57fa999abe8dc557003482444df2a189) to disable the scaling, which otherwise aims to compensate volume loss during smoothing.

-   Deprecated the overloads of functions [`CGAL::Polygon_mesh_processing::triangulate_hole()`](https://doc.cgal.org/5.6/Polygon_mesh_processing/group__PMP__hole__filling__grp.html#ga3abdf2d0558822e85f060966b69cae98), [`CGAL::Polygon_mesh_processing::triangulate_and_refine_hole()`](https://doc.cgal.org/5.6/Polygon_mesh_processing/group__PMP__hole__filling__grp.html#ga9868fac4d9dca77462ad7828bc99d8a1), and [`CGAL::Polygon_mesh_processing::triangulate_refine_and_fair_hole()`](https://doc.cgal.org/5.6/Polygon_mesh_processing/group__PMP__hole__filling__grp.html#ga18eac756a8f8e5d5f73e645fd4e26cad) which have output iterators for vertices and faces as parameter. They are replaced by overloads with two additional named parameters.

### [2D Convex Hulls](https://doc.cgal.org/5.6/Manual/packages.html#PkgConvexHull2)

-   **Breaking change**: The concept [`ConvexHullTraits_2`](https://doc.cgal.org/5.6/Convex_hull_2/classConvexHullTraits__2.html) no longer requires the functor `Less_signed_distance_to_line_2`, but requires the functor `Compare_signed_distance_to_line_2` instead.
-   The long-deprecated classes `Convex_hull_projective_xy_traits_2`, `Convex_hull_projective_xz_traits_2`, and `Convex_hull_projective_yz_traits_2` have been removed. Users should use [`Projection_traits_xy_3`](https://doc.cgal.org/5.6/Kernel_23/classCGAL_1_1Projection__traits__xy__3.html), [`Projection_traits_xz_3`](https://doc.cgal.org/5.6/Kernel_23/classCGAL_1_1Projection__traits__xz__3.html), and [`Projection_traits_yz_3`](https://doc.cgal.org/5.6/Kernel_23/classCGAL_1_1Projection__traits__yz__3.html) instead.

### [2D Triangulations](https://doc.cgal.org/5.6/Manual/packages.html#PkgTriangulation2)

-   Added the function [`CGAL::mark_domain_in_triangulation()`](https://doc.cgal.org/5.6/Triangulation_2/group__PkgTriangulation2Miscellaneous.html#ga0409755d0eb89100810230443a85e7eb) to mark faces connected with non-constrained edges as inside of the domain based on the nesting level.

### [2D Conforming Triangulations and Meshes](https://doc.cgal.org/5.6/Manual/packages.html#PkgMesh2)

-   Added new overloads to the function [`write_VTU()`](https://doc.cgal.org/5.6/Mesh_2/group__PkgMesh2IO.html), with property maps for specifying the domain.
-   Deprecated usage of boost parameters in favor of function named parameters in [`CGAL::lloyd_optimize_mesh_2()`](https://doc.cgal.org/5.6/Mesh_2/group__PkgMesh2Functions.html#gafeaf59d3fa014da287f8514913b38d05).
-   Deprecated two overloads of the function [`refine_Delaunay_mesh()`](https://doc.cgal.org/5.6/Mesh_2/group__PkgMesh2Functions.html), and replaced them with versions using function named parameters.

### [2D Hyperbolic Triangulations](https://doc.cgal.org/5.6/Manual/packages.html#PkgHyperbolicTriangulation2)

-   **Breaking change**: the concept [`HyperbolicTriangulationFaceBase_2`](https://doc.cgal.org/5.6/Hyperbolic_triangulation_2/classHyperbolicTriangulationFaceBase__2.html) has been modified to better reflect the triangulation's requirements and avoid a conflict with the requirements described by the concept `TriangulationDataStructure_2::Face`. The model [`CGAL::Hyperbolic_triangulation_face_base_2`](https://doc.cgal.org/5.6/Hyperbolic_triangulation_2/classCGAL_1_1Hyperbolic__triangulation__face__base__2.html) has been adapted correspondingly.

### [3D Simplicial Mesh Data Structure](https://doc.cgal.org/5.6/Manual/packages.html#PkgSMDS3) (new package)

-   This new package wraps all the existing code that deals with a [`MeshComplex_3InTriangulation_3`](https://doc.cgal.org/5.6/SMDS_3/classMeshComplex__3InTriangulation__3.html) to describe 3D simplicial meshes, and makes the data structure independent from the [tetrahedral mesh generation](https://doc.cgal.org/5.6/Manual/packages.html#PkgMesh3) package.

### [3D Mesh Generation](https://doc.cgal.org/5.6/Manual/packages.html#PkgMesh3)

-   Added two new named parameters to the named constructor [`CGAL::create_labeled_image_mesh_domain()`](https://doc.cgal.org/5.6/Mesh_3/classCGAL_1_1Labeled__mesh__domain__3.html#aec3f58e9883a8036a1b3e379df7d8fa9) for automatic detection and protection of 1D-curves that lie at the intersection of three or more subdomains extracted from labeled images.
-   Added [`CGAL::Sizing_field_with_aabb_tree`](https://doc.cgal.org/5.6/Mesh_3/structCGAL_1_1Sizing__field__with__aabb__tree.html), a geometry-aware sizing field for feature edges in polyhedral domains.
-   Added new meshing criterion [`edge_min_size`](https://doc.cgal.org/5.6/Mesh_3/classCGAL_1_1Mesh__criteria__3.html#a5f1c2649cb7ea346a3b6a2a8724b4df1) to avoid subdividing sharp edges that are shorter than a prescribed size bound.
-   Added new meshing criteria [`facet_min_size`](https://doc.cgal.org/5.6/Mesh_3/classCGAL_1_1Mesh__criteria__3.html#a5f1c2649cb7ea346a3b6a2a8724b4df1) and [`cell_min_size`](https://doc.cgal.org/5.6/Mesh_3/classCGAL_1_1Mesh__criteria__3.html#a5f1c2649cb7ea346a3b6a2a8724b4df1) to prevent Delaunay refinement from creating simplices smaller than a prescribed bound.
-   Deprecated usage of boost parameters in favor of function named parameters.

### [3D Periodic Mesh Generation](https://doc.cgal.org/5.6/Manual/packages.html#PkgPeriodic3Mesh3)

-   Periodic Mesh Generation now supports non-cubic domains.
-   Deprecated usage of boost parameters in favor of function named parameters.

### [Surface Mesh Simplification](https://doc.cgal.org/5.6/Manual/packages.html#PkgSurfaceMeshSimplification)
-   The stop predicates [`Count_stop_predicate`](https://doc.cgal.org/5.6/Surface_mesh_simplification/classCGAL_1_1Surface__mesh__simplification_1_1Count__stop__predicate.html) and [`Count_ratio_stop_predicate`](https://doc.cgal.org/5.6/Surface_mesh_simplification/classCGAL_1_1Surface__mesh__simplification_1_1Count__ratio__stop__predicate.html) are renamed to [`Edge_count_stop_predicate`](https://doc.cgal.org/5.6/Surface_mesh_simplification/classCGAL_1_1Surface__mesh__simplification_1_1Edge__count__stop__predicate.html) and [`Edge_count_ratio_stop_predicate`](https://doc.cgal.org/5.6/Surface_mesh_simplification/classCGAL_1_1Surface__mesh__simplification_1_1Edge__count__ratio__stop__predicate.html). Older versions have been deprecated.
-   Introduced [`Face_count_stop_predicate`](https://doc.cgal.org/5.6/Surface_mesh_simplification/classCGAL_1_1Surface__mesh__simplification_1_1Face__count__stop__predicate.html) and [`Face_count_ratio_stop_predicate`](https://doc.cgal.org/5.6/Surface_mesh_simplification/classCGAL_1_1Surface__mesh__simplification_1_1Face__count__ratio__stop__predicate.html), which can be used to stop the simplification algorithm based on a desired number of faces in the output, or a ratio between input and output face numbers.

### [2D Regularized Boolean Set Operations](https://doc.cgal.org/5.6/Manual/packages.html#PkgBooleanSetOperations2)
-   Exposed all required member functions of the [`GeneralPolygonWithHoles_2`](https://doc.cgal.org/5.6/Polygon/classGeneralPolygonWithHoles__2.html) concept (e.g., [`clear_outer_boundary()`](https://doc.cgal.org/5.6/Polygon/classGeneralPolygonWithHoles__2.html#a9f5f035047505a2ccab3e68770f51bc6), [`clear_holes()`](https://cgal.geometryfactory.com/CGAL/doc/master/Polygon/classGeneralPolygonWithHoles__2.html#a2a507be648f127ac605da8c670ea2580), and [`clear()`](https://doc.cgal.org/5.6/Polygon/classGeneralPolygonWithHoles__2.html#a2ca4d9b43cc9216c1b2cdb080a915944) ).
