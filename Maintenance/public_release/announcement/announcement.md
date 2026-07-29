The CGAL Open Source Project is pleased to announce the release 6.2 of CGAL, the Computational Geometry Algorithms Library.

Besides fixes and general enhancement to existing packages, the following has changed since [CGAL 6.1]:

### General Changes

- The new list of supported compilers is:
  - Visual C++ 15.9, 16.10, 17.14, 18.0 (from Visual Studio 2017, 2019, 2022, and 2026) or later
  - Gnu g++ 13.3.0 or later (on Linux)
  - LLVM Clang version 22.0.0 or later (on Linux)
  - Apple Clang compiler versions 14.0.0 or later (on macOS)
- The minimal supported version of Boost is still 1.74.0, but only versions 1.79.0 and later were tested.

### [2D Alpha Wrapping](https://doc.cgal.org/6.2/Manual/packages.html#PkgAlphaWrap2) (new package)

- This component takes a polygon soup, a 2D segment soup, and/or a 2D point set as input, and generates a valid (watertight, intersection-free and 1-manifold) multi-polygon that strictly encloses the input. The algorithm proceeds by shrink-wrapping and refining a 2D Delaunay triangulation starting from a loose bounding box of the input. Two user-defined parameters, alpha and offset, offer control over the maximum size of cavities where the shrink-wrapping process can enter, and the tightness of the final polygon(s) to the input, respectively. Once combined, these parameters provide a means to trade fidelity to the input for complexity of the output.

  See also the associated [news entry](https://www.cgal.org/2026/03/30/alpha_wrap_2/).

### [Generalized Barycentric Coordinates 3](https://doc.cgal.org/6.2/Manual/packages.html#PkgBarycentricCoordinates3) (new package)

- This package provides functions to compute various types of generalized barycentric coordinates (Wachspress, mean value, discrete harmonic and tetrahedron coordinates) for points located inside closed convex 3D polyhedra.

### [Polygon Mesh Processing](https://doc.cgal.org/6.2/Manual/packages.html#PkgPolygonMeshProcessing) (major changes)

- The “Polygon Mesh Processing” package has been reorganized into several packages. “Polygon Mesh Processing” retains the core functionalities, while advanced and specialized features have been moved to dedicated packages:
  - [Boolean Operations On Meshes](https://doc.cgal.org/6.2/Manual/packages.html#PkgPMPBooleanOperations): algorithms for Boolean operations on polygon meshes; clipping, splitting, and slicing with planes, boxes, or other meshes; and kernel computations.
  - [Meshing and Remeshing of Polygon Meshes](https://doc.cgal.org/6.2/Manual/packages.html#PkgPMPRemeshing): algorithms for meshing and remeshing, such as triangulation, refinement, simplification, optimization, and smoothing.
  - [Polygon Mesh Repair](https://doc.cgal.org/6.2/Manual/packages.html#PkgPMPMeshRepair): tools for detecting and correcting combinatorial and geometric defects in polygon meshes and polygon soups, including face orientation, hole filling, removal of degeneracies, and boundary stitching. This split does not induce any breaking change and is fully transparent: header includes such as `#include <CGAL/Polygon_mesh_processing/XXX.h>` do not need to be changed and will include the appropriate header from the new packages.

### [2D and 3D Linear Geometry Kernel](https://doc.cgal.org/6.2/Manual/packages.html#PkgKernel23)

- **Breaking change**: The functors [`Do_intersect_2`](https://doc.cgal.org/6.2/Kernel_23/classKernel.html#a66956c9eb7a1e2042326b051ecef8020) and [`Do_intersect_3`](https://doc.cgal.org/6.2/Kernel_23/classKernel.html#adf4b45486ddac1219aa419fc62c38e1b) for intersection now correctly honor the definition of a `Circle_2` (resp. `Sphere_3`) as a 1-manifold (resp. 2-manifold) without borders. In other words, a circle is not a disk, and a sphere is not a ball. Consequently, `Circle_2`/`Segment_2`, `Sphere_3`/`Bbox_3`, `Sphere_3`/`Iso_cuboid_3` no longer consider inclusion as an intersection. This behavior is consistent with other intersection combinations involving `Circle_2` and `Sphere_3`. The former behavior can be reproduced with:
  - [`!Has_on_unbounded_side_2::operator(Circle_2, Segment_2)`](https://doc.cgal.org/6.2/Kernel_23/classKernel.html#ab041aafd45ce93630bfe77cec24e24f3),
  - [`!Has_on_unbounded_side_2::operator(Circle_2, Iso_rectangle_2)`](https://doc.cgal.org/6.2/Kernel_23/classKernel.html#ab041aafd45ce93630bfe77cec24e24f3), or
  - [`!Has_on_unbounded_side_3::operator(Sphere_3, Iso_cuboid_3)`](https://doc.cgal.org/6.2/Kernel_23/classKernel.html#ada27f27602757133ea7b9814ada00ba3)
- Added a new concept, [`CompareProjectionAlongDirection_3`](https://doc.cgal.org/6.2/Kernel_23/classKernel_1_1CompareProjectionAlongDirection__3.html), to the [`Kernel`](https://doc.cgal.org/6.2/Kernel_23/classKernel.html) concept to compare the order of point projections onto a line. Corresponding functors in the model, [`Compare_projection_along_direction_3`](https://doc.cgal.org/6.2/Kernel_23/classKernel.html#a4b7686814790650efe8f8b6612af3b80), and free function, [`compare_projection_along_direction()`](https://doc.cgal.org/6.2/Kernel_23/group__kernel__global__function.html#ga2ba880adb55cea48d3568ccd1f4765a4), have been added.

### [Intersecting Sequences of dD Iso-oriented Boxes](https://doc.cgal.org/6.2/Manual/packages.html#PkgBoxIntersectionD)

- The function [`CGAL::box_intersection_d()`](https://doc.cgal.org/6.2/Box_intersection_d/group__PkgBoxIntersectionD__box__intersection__d.html) now supports ranges with different box types.
- Added the box type [`CGAL::Box_with_info_d<FT, dim, Info>`](https://doc.cgal.org/6.2/Box_intersection_d/classCGAL_1_1Box__intersection__d_1_1Box__with__info__d.html), which stores a variable of type `Info` accessible via the class member function [`info()`](https://doc.cgal.org/6.2/Box_intersection_d/classCGAL_1_1Box__intersection__d_1_1Box__with__info__d.html#a5858777d3531cfc13bdc0fea067a7532).

### [2D Arrangements](https://doc.cgal.org/6.2/Manual/packages.html#PkgArrangementOnSurface2)

- Introduced a Geometry Traits concept for arrangement on surfaces that enables the provision of the disconnected portions of an approximation of a curve within a given bounding box.
- The class [`CGAL::Arr_linear_traits_2`](https://doc.cgal.org/6.2/Arrangement_on_surface_2/classCGAL_1_1Arr__linear__traits__2.html) is now a model of the new concept.
- Added overloads of [`CGAL::draw(Arrangement_on_surface_2& arr, Bbox& bbox, ...)`](https://doc.cgal.org/6.2/Arrangement_on_surface_2/group__PkgArrangementOnSurface2Draw.html) that enable the drawing of arrangements induced by unbounded curves.
- Introduced the concept `AosTraits::Do_intersect_2`. A model of this concept must provide an operator that accepts two `x`-monotone curves and a boolean flag that indicates whether common endpoints should be considered or ignored. The operator determines whether the curves intersect, and it can be used with an inexact-construction kernel.

### [2D Intersection of Curves](https://doc.cgal.org/6.2/Manual/packages.html#PkgSurfaceSweep2)

- The function [`CGAL::do_curves_intersect()`](https://doc.cgal.org/6.2/Surface_sweep_2/group__PkgSurfaceSweep2Ref.html#gadb96b95f091371e834e3f7b4d24f7bb4), which assumed open curves, has been deprecated and replaced by the function [`CGAL::Surface_sweep_2::do_intersect()`](https://doc.cgal.org/6.2/Surface_sweep_2/group__PkgSurfaceSweep2Ref.html#ga6de86b9ce9ed5537d371142716cf8adf). Notice: (i) the introduction of the new namespace `Surface_sweep_2`, and (ii) the addition of the `closed` parameter, which defaults to `true`. To match the behavior of the deprecated function, set `closed` to false.

### [2D Regularized Boolean Set-Operations](https://doc.cgal.org/6.2/Manual/packages.html#PkgBooleanSetOperations2)

- Optimized [`do_intersect()`](https://doc.cgal.org/6.2/Boolean_set_operations_2/group__boolean__do__intersect.html):
  - (i)  made it robust even with an inexact-predicates kernel, and
  - (ii)  made it quit once an intersection is detected. (In the past, the intersection was computed in one phase and examined in a subsequent phase.) This optimization somehow breaks backward compatibility as follows: the variants of the free function `do_intersect()` that accept a third optional parameter, namely `UsePolylines`, which determines whether the boundaries of the input polygons are treated as cyclic sequences of (`x`-monotone) segments or as a cyclic sequences of (`x`-monotone) polylines, do not accept this third parameter any longer. (This third optional parameter was introduced a few years ago, and now abandoned only for `do_intersect()`.)

### [3D Convex hull](https://doc.cgal.org/6.2/Manual/packages.html#PkgConvexHull3)

- Added the function [`CGAL::Convex_hull_3::do_intersect()`](https://doc.cgal.org/6.2/Convex_hull_3/group__PkgConvexHull3Intersections.html#gaf41281925d9bf8b5b0c6e13851d900c7), which can be used to test for intersection between two convex hulls.
- Added the function [`CGAL::extreme_point_3()`](https://doc.cgal.org/6.2/Convex_hull_3/group__PkgConvexHull3Queries.html#ga329e2cdfdc3bdf89a47073c277ce202d), which returns the farthest point of a convex hull in a given direction.
- Added the class [`CGAL::Convex_hull_hierarchy_3`](https://doc.cgal.org/6.2/Convex_hull_3/structCGAL_1_1Convex__hull__hierarchy__3.html), which represents a convex hull, and is optimized for the functions above.

### [2D Triangulations](https://doc.cgal.org/6.2/Manual/packages.html#PkgTriangulation2)

- Added the function [`insert_unique_constraints()`](https://doc.cgal.org/6.2/Triangulation_2/classCGAL_1_1Constrained__Delaunay__triangulation__2.html#a228faa9b87674183d88e3858851e9747) to the class [`CGAL::Constrained_Delaunay_triangulation_2`](https://doc.cgal.org/6.2/Triangulation_2/classCGAL_1_1Constrained__Delaunay__triangulation__2.html) identical to the function [`insert_constraints()`](https://doc.cgal.org/6.2/Triangulation_2/classCGAL_1_1Constrained__Delaunay__triangulation__2.html#a2fa90339d0c18d07c60857fe42acab2a) except that it removes duplicated constraints before inserting them in the triangulation.

### [2D Conforming Triangulations and Meshes](https://doc.cgal.org/6.2/Manual/packages.html#PkgMesh2)

- The implementation is now more robust to almost degenerate inputs, such as polygons with microscopic edges or almost collinear points.
- **Breaking change**: The concept [`DelaunayMeshTraits_2`](https://doc.cgal.org/6.2/Mesh_2/classDelaunayMeshTraits__2.html) now requires the functor `Construct_bbox_2`.

### [Tetrahedral Mesh Generation](https://doc.cgal.org/6.2/Manual/packages.html#PkgMesh3)

- **Breaking change**: Removed the class template `CGAL::Implicit_vector_to_labeling_function_wrapper`, as well as the constructor of [`CGAL::Polyhedral_mesh_domain_with_features_3`](https://doc.cgal.org/6.2/Mesh_3/classCGAL_1_1Polyhedral__mesh__domain__with__features__3.html) that takes a filename as parameter, which were deprecated since CGAL-4.5.
- **Breaking change**: Added the requirement for a nested type `Iso_cuboid_3` to the concept [`BisectionGeometricTraits_3`](https://doc.cgal.org/6.2/Mesh_3/classBisectionGeometricTraits__3.html).
- Protection of sharp edges (also known as “feature line”) is now significantly faster.

### [dD Triangulations](https://doc.cgal.org/6.2/Manual/packages.html#PkgTriangulations)

- Computation of convex hulls in high dimensions is now significantly faster.

### [Convex Decomposition of Polyhedra](https://doc.cgal.org/6.2/Manual/packages.html#PkgConvexDecomposition3)

- Added the function [`CGAL::approximate_convex_decomposition()`](https://doc.cgal.org/6.2/Convex_decomposition_3/group__PkgConvexDecomposition3Ref.html#ga98ed6393d56fc2eac2ad82a55af84b2f), which computes a set of convex volumes that cover an input mesh.

### [Polygon Mesh Processing](https://doc.cgal.org/6.2/Manual/packages.html#PkgPolygonMeshProcessing)

- **Breaking change**: The header `CGAL/Polygon_mesh_processing/border.h` has been deprecated and its content ([`CGAL::border_halfedges()`](https://doc.cgal.org/6.2/BGL/group__PkgBGLHelperFct.html#gad9e4e6644bc6e2c79101b25e280ea633) and [`CGAL::extract_boundary_cycles()`](https://doc.cgal.org/6.2/BGL/group__PkgBGLHelperFct.html#ga28bba3efbd0352497f28c3c8185be09f)) have been moved to a new header: `CGAL/boost/graph/border.h`.
- Vertex normal computation is now significantly faster for high-degree vertices.

### [Polygon Mesh Processing (Boolean Operations on Meshes)](https://doc.cgal.org/6.2/Manual/packages.html#PkgPMPBooleanOperations)

- Added the function [`CGAL::Polygon_mesh_processing::kernel()`](https://doc.cgal.org/6.2/PMP_Boolean_operations/group__PMP__kernel__grp.html#ga778059f8efd21fb9a9ffc7ede16c5abc), which can be used to compute the kernel of a polygon mesh, and is significantly faster than halfspace intersections.
- Added the function [`CGAL::Polygon_mesh_processing::has_empty_kernel()`](https://doc.cgal.org/6.2/PMP_Boolean_operations/group__PMP__kernel__grp.html#gaf5b25f2dc1e3d96e9f81ed182f9600de), which can be used to determine if the kernel of a polygon mesh is empty.
- Added the function [`CGAL::Polygon_mesh_processing::kernel_point()`](https://doc.cgal.org/6.2/PMP_Boolean_operations/group__PMP__kernel__grp.html#gaea657a1fbac5d35c3d51735c3ec2fede), which can be used to compute a single point inside the kernel of a polygon mesh.
- Added an optional parameter `use_convex_specialization` to the functions [`CGAL::Polygon_mesh_processing::clip()`](https://doc.cgal.org/6.2/PMP_Boolean_operations/group__PMP__clip__grp.html#ga2c73d3460872e601f84a03f58dd069ae) and [`CGAL::Polygon_mesh_processing::refine_with_plane()`](https://doc.cgal.org/6.2/PMP_Boolean_operations/group__PMP__slicing__grp.html#gacb9d68fa4dea8fd03ec53b56a16d6fc6), which can be used to specify that a given input is convex, resulting in much faster computation.
- The function [`CGAL::Polygon_mesh_processing::surface_intersection()`](https://doc.cgal.org/6.2/PMP_Boolean_operations/group__PMP__corefinement__grp.html#ga6e6c4a724ce19e7a207de56f3a7408ab), which can be used to compute the polylines that represent the intersection between two polyhedral surfaces has been renamed to [`CGAL::Polygon_mesh_processing::intersection_polylines()`](https://doc.cgal.org/6.2/PMP_Boolean_operations/group__PMP__corefinement__grp.html#ga7da80d8488d7a9139dcabceeb601072f).

### [Polygon Mesh Processing (Meshing and Remeshing of Polygon Meshes)](https://doc.cgal.org/6.2/Manual/packages.html#PkgPMPRemeshing)

- **Breaking change**: Update the visitor concepts [`PMPTriangulateFaceVisitor`](https://doc.cgal.org/6.2/Polygon_mesh_processing/classPMPTriangulateFaceVisitor.html) and [`PMPHolefillingVisitor`](https://doc.cgal.org/6.2/Polygon_mesh_processing/classPMPHolefillingVisitor.html) to request functions [`accept_face()`](https://doc.cgal.org/6.2/Polygon_mesh_processing/classPMPTriangulateFaceVisitor.html#a62ee542bc34faafb193147c746e07a69) and [`accept_triangle()`](https://doc.cgal.org/6.2/Polygon_mesh_processing/classPMPHolefillingVisitor.html#a273b240262220abe572e42d8daef0d7a), respectively. These functions can be used to tweak the way faces and holes are triangles by black listing some candidate triangles. User visitors inheriting from the default visitors do not require any update.

### [Polygon Mesh Processing (Mesh Repair)](https://doc.cgal.org/6.2/Manual/packages.html#PkgPMPMeshRepair)

- Added the named parameter `erase_policy` to [`CGAL::Polygon_mesh_processing::repair_polygon_soup()`](https://doc.cgal.org/6.2/PMP_Mesh_repair/group__PMP__combinatorial__repair__grp.html#ga3b35133783759402828325b91ab559cc) and [`CGAL::Polygon_mesh_processing::merge_duplicate_polygons_in_polygon_soup()`](https://doc.cgal.org/6.2/PMP_Mesh_repair/group__PMP__combinatorial__repair__grp.html#ga1f215926ed8794db827e2993d2960870). This parameter offers three policies: (i) erase all duplicates polygons,
  2)  keep one of the duplicates (default), and (iii) keep one iff the number of duplicates is odd. The named parameter `erase_all_duplicates` is now deprecated.

### [Quadtrees, Octrees, and Orthtrees](https://doc.cgal.org/6.2/Manual/packages.html#PkgOrthtree)

- Added the member function [`intersected_nodes()`](https://doc.cgal.org/6.2/Orthtree/classCGAL_1_1Orthtree.html#a595f6d7c735f4c60c2c65e5ed2ef293c) [with an intersection functor](https://doc.cgal.org/6.2/Orthtree/classCGAL_1_1Orthtree.html#a18a2ac4bb94e3dbc17e556e833d2fe68) and [a convenience overload for a ball query](https://doc.cgal.org/6.2/Orthtree/classCGAL_1_1Orthtree.html#af685c4789f52eac1946add391b849777) to the class [`CGAL::Orthree`](https://doc.cgal.org/6.2/Orthtree/classCGAL_1_1Orthtree.html).

### [Linear Cell Complex](https://doc.cgal.org/6.2/Manual/packages.html#PkgLinearCellComplex)

- The following import functions have been deprecated and renamed for better naming clarity and consistency:
  - [`import_from_plane_graph()`](https://doc.cgal.org/6.2/Linear_cell_complex/group__PkgLinearCellComplexConstructions.html#ga2e694849af5a899bb8545ff7ea8ab2a8) → [`read_plane_graph_in_lcc()`](https://doc.cgal.org/6.2/Linear_cell_complex/group__PkgLinearCellComplexConstructions.html#ga909e9e304b69e21ae3cddc0f5bb7b94a)
  - `import_from_polyhedron_3()` → [`polyhedron_3_to_lcc()`](https://doc.cgal.org/6.2/Linear_cell_complex/group__PkgLinearCellComplexConstructions.html#gac790497ffcb75b93dc6994d8239628a9)
  - `import_from_triangulation_3()` → [`triangulation_3_to_lcc()`](https://doc.cgal.org/6.2/Linear_cell_complex/group__PkgLinearCellComplexConstructions.html#gae0d938bab5a5da62716b67a69f4cd520) The old function names are still available but marked as deprecated for backward compatibility.
- Two functions have been added to read and write [3D linear cell complexes](https://doc.cgal.org/6.2/Linear_cell_complex/classLinearCellComplex.html) (dimension 3, ambient dimension 3) using VTK file formats (legacy ASCII):
  - [`CGAL::IO::read_VTK<LCC>()`](https://doc.cgal.org/6.2/Linear_cell_complex/group__PkgLinearCellComplexRefIOVTK.html#ga5c283bb942b300eb47c956ae63157507), and
  - [`CGAL::IO::write_VTK<LCC>()`](https://doc.cgal.org/6.2/Linear_cell_complex/group__PkgLinearCellComplexRefIOVTK.html#ga857ed2854ab8fe2c65126fb066dd5e07) These functions support per-vertex and per-volume scalar fields and handle various VTK cell types.

### [Shape Detection](https://doc.cgal.org/6.2/Manual/packages.html#PkgShapeDetection)

- Added the region type [`CGAL::Shape_detection::Polygon_mesh::Plane_face_region`](https://doc.cgal.org/6.2/Shape_detection/classCGAL_1_1Shape__detection_1_1Polygon__mesh_1_1Plane__face__region.html) that extends the support plane of the seed face without refitting the plane to the region
- Added the region type [`CGAL::Shape_detection::Polygon_mesh::Line_segment_region`](https://doc.cgal.org/6.2/Shape_detection/classCGAL_1_1Shape__detection_1_1Segment__set_1_1Line__segment__region.html) that extends the support line of the seed segment without refitting the line to the region
- Added the sorting [`CGAL::Shape_detection::Polygon_mesh::Face_area_sorting`](https://doc.cgal.org/6.2/Shape_detection/classCGAL_1_1Shape__detection_1_1Polygon__mesh_1_1Face__area__sorting.html) that provides a sorting of faces in descending order of area
- Added the sorting [`CGAL::Shape_detection::Polygon_mesh::Segment_length_sorting`](https://doc.cgal.org/6.2/Shape_detection/classCGAL_1_1Shape__detection_1_1Segment__set_1_1Segment__length__sorting.html) that provides a sorting of segments in descending order of length
- Added the optional named parameter `face_normal_map` to [`CGAL::Shape_detection::Polygon_mesh::Least_squares_plane_fit_region`](https://doc.cgal.org/6.2/Shape_detection/classCGAL_1_1Shape__detection_1_1Polygon__mesh_1_1Least__squares__plane__fit__region.html), [`CGAL::Shape_detection::Polygon_mesh::Plane_face_region`](https://doc.cgal.org/6.2/Shape_detection/classCGAL_1_1Shape__detection_1_1Polygon__mesh_1_1Plane__face__region.html); and [`CGAL::Polygon_mesh_processing::region_growing_of_planes_on_faces()`](https://doc.cgal.org/6.2/Polygon_mesh_processing/group__PMP__detect__features__grp.html#ga50dcd2f6295f584d2e378b57290ae2af) to allow the use of face normal property maps instead of calculating the face normals.
- The [`CGAL::Shape_detection::Polygon_mesh::Plane_face_region`](https://doc.cgal.org/6.2/Shape_detection/classCGAL_1_1Shape__detection_1_1Polygon__mesh_1_1Plane__face__region.html) and [`CGAL::Shape_detection::Polygon_mesh::Face_area_sorting`](https://doc.cgal.org/6.2/Shape_detection/classCGAL_1_1Shape__detection_1_1Polygon__mesh_1_1Face__area__sorting.html) are now used as the default in [`CGAL::Polygon_mesh_processing::region_growing_of_planes_on_faces()`](https://doc.cgal.org/6.2/Polygon_mesh_processing/group__PMP__detect__features__grp.html#ga50dcd2f6295f584d2e378b57290ae2af)

### [Surface Mesh Simplification](https://doc.cgal.org/6.2/Manual/packages.html#PkgSurfaceMeshSimplification)

- Added the class [`CGAL::Surface_mesh_simplification::GarlandHeckbert_plane_and_line_policies`](https://doc.cgal.org/6.2/Surface_mesh_simplification/classCGAL_1_1Surface__mesh__simplification_1_1GarlandHeckbert__plane__and__line__policies.html), which improves other Garland-Heckbert strategies near sharp features.
- **Breaking change**: [`CGAL::Surface_mesh_simplification::GarlandHeckbert_policies.h`](https://doc.cgal.org/6.2/Surface_mesh_simplification/group__PkgSurfaceMeshSimplificationRef.html#ga3f012cf5dc52d6e5ac2aea1106e4cad7) is is no longer deprecated and has become an alias for the current best Garland-Heckbert policy, which is currently `CGAL::Surface_mesh_simplification::GarlandHeckbert_plane_and_line_policies`.

### [Geometric Object Generators](https://doc.cgal.org/6.2/Generator/index.html)

- Added three new point generators:
  - [`CGAL::Random_points_on_segment_3`](https://doc.cgal.org/6.2/Generator/classCGAL_1_1Random__points__on__segment__3.html),
  - [`CGAL::Random_points_in_triangle_soup_3`](https://doc.cgal.org/6.2/Generator/structCGAL_1_1Random__points__in__triangle__soup__3.html), and
  - [`CGAL::Random_points_on_graph_edges_3`](https://doc.cgal.org/6.2/Generator/structCGAL_1_1Random__points__on__graph__edges__3.html)
- Added the function `point_and_support()` to several generators to get the last point generated point and the support element containing it. Affected generators are:
  - [`CGAL::Random_points_in_triangles_2`](https://doc.cgal.org/6.2/Generator/structCGAL_1_1Random__points__in__triangles__2.html),
  - [`CGAL::Random_points_in_triangles_3`](https://doc.cgal.org/6.2/Generator/structCGAL_1_1Random__points__in__triangles__3.html),
  - [`CGAL::Random_points_in_triangle_mesh_3`](https://doc.cgal.org/6.2/Generator/structCGAL_1_1Random__points__in__triangle__mesh__3.html),
  - [`CGAL::Random_points_in_tetrahedral_mesh_boundary_3`](https://doc.cgal.org/6.2/Generator/classCGAL_1_1Random__points__in__tetrahedral__mesh__boundary__3.html),
  - [`CGAL::Random_points_in_tetrahedral_mesh_3`](https://doc.cgal.org/6.2/Generator/classCGAL_1_1Random__points__in__tetrahedral__mesh__3.html),
  - [`CGAL::Random_points_in_triangle_soup_3`](https://doc.cgal.org/6.2/Generator/structCGAL_1_1Random__points__in__triangle__soup__3.html),
  - and [`CGAL::Random_points_on_graph_edges_3`](https://doc.cgal.org/6.2/Generator/structCGAL_1_1Random__points__on__graph__edges__3.html).

### [STL Extension](https://doc.cgal.org/6.2/Manual/packages.html#PkgSTLExtension)

- Added new debugging utility for finding minimal failing test cases:
  - [`CGAL::bisect_failures`](https://doc.cgal.org/6.2/STL_Extension/group__PkgSTLExtensionUtilities.html#gad56d708faac9063d3de057394a6d8c52) is a template function that uses bisection to isolate minimal failing subsets from complex input data.

### [Stream Support](https://doc.cgal.org/6.2/Manual/packages.html#PkgStreamSupport)

- Added new stream formatting capabilities for improved debugging and logging:
  - [`CGAL::IO::Basic_indenting_streambuf`](https://doc.cgal.org/6.2/Stream_support/classCGAL_1_1IO_1_1Basic__indenting__streambuf.html) and [`CGAL::IO::Basic_indenting_stream_guard`](https://doc.cgal.org/6.2/Stream_support/classCGAL_1_1IO_1_1Basic__indenting__stream__guard.html) for automatic indentation of output streams.
  - [`CGAL::IO::Basic_color_streambuf`](https://doc.cgal.org/6.2/Stream_support/classCGAL_1_1IO_1_1Basic__color__streambuf.html) and [`CGAL::IO::Basic_color_stream_guard`](https://doc.cgal.org/6.2/Stream_support/classCGAL_1_1IO_1_1Basic__color__stream__guard.html) for ANSI color support in terminal output.

[CGAL 6.1]: https://github.com/CGAL/cgal/releases/tag/v6.1
