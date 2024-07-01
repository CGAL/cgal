The CGAL Open Source Project is pleased to announce the release 6.0 Beta 1 of CGAL, the Computational Geometry Algorithms Library.

CGAL version 6.0 Beta 1 is a public testing release. It should provide a solid ground to report bugs that need to be tackled before the release of the final version of CGAL 6.0 in July 2024.

Besides fixes and general enhancement to existing packages, the following has changed since CGAL 5.6:

### General Changes

- CGAL 6.0 is the first release of CGAL that requires a C++ compiler with the support of C++17 or later. The new list of supported compilers is:
  - Visual C++ 15.9, 16.10, 17.0 (from Visual Studio 2017, 2019 and 2022) or later
  - Gnu g++ 11.4.0 or later (on Linux or macOS)
  - LLVM Clang version 15.0.7 or later (on Linux)
  - Apple Clang compiler versions 10.0.1, 12.0.5, and 15.0.0 (on macOS)
- The minimal supported version of Boost is now 1.72.0.
- The CGAL `Core` library is no longer based on GMP, but on [Boost.Multiprecision](https://www.boost.org/doc/libs/1_72_0/libs/multiprecision/doc/html/index.html). Either GMP backend or Boost backend can be used.
- All demos are now based on Qt6.
- **Breaking change**: The CMake file `UseCGAL.cmake` has been removed from CGAL. Usages of the CMake variables `${CGAL_USE_FILE}` and `${CGAL_LIBRARIES}` must be replaced by a link to the imported target `CGAL::CGAL`, for example: `target_link_library(your_target PRIVATE CGAL::CGAL)`.

### [Kinetic Space Partition](https://doc.cgal.org/6.0/Manual/packages.html#PkgKineticSpacePartition) (new package)

- This package implements kinetic space partition: based on a set of planar input shapes, the bounding box of the input data is split into convex volumes. The complexity of the partition can be adjusted with a single parameter.

### [Kinetic Surface Reconstruction](https://doc.cgal.org/6.0/Manual/packages.html#PkgKineticSurfaceReconstruction) (new package)

- The package implements a piece-wise planar surface reconstruction pipeline from point clouds combining methods from the [Shape Detection](https://doc.cgal.org/6.0/Manual/packages.html#PkgShapeDetection), [Shape Regularization](https://doc.cgal.org/6.0/Manual/packages.html#PkgShapeRegularization) and [Kinetic Shape Partition](https://doc.cgal.org/6.0/Manual/packages.html#PkgKineticSpacePartition) packages and graph-cut to reconstruct surfaces from point clouds.

### [Basic Viewer](https://doc.cgal.org/6.0/Basic_viewer/index.html#Chapter_Basic_viewer) (new package)

- The basic viewer package provides interactive visualization for most CGAL packages, such as [2D Arrangements](https://doc.cgal.org/6.0/Manual/packages.html#PkgArrangementOnSurface2), [2D Regularized Boolean Set-Operations](https://doc.cgal.org/6.0/Manual/packages.html#PkgBooleanSetOperations2), [Linear Cell Complex](https://doc.cgal.org/6.0/Manual/packages.html#PkgLinearCellComplex), [3D Boolean Operations on Nef Polyhedra](https://doc.cgal.org/6.0/Manual/packages.html#PkgNef3), [2D Periodic Triangulations](https://doc.cgal.org/6.0/Manual/packages.html#PkgPeriodic2Triangulation2), [3D Point Set](https://doc.cgal.org/6.0/Manual/packages.html#PkgPointSet3), [2D Polygons](https://doc.cgal.org/6.0/Manual/packages.html#PkgPolygon2), [3D Polyhedral Surface](https://doc.cgal.org/6.0/Manual/packages.html#PkgPolyhedron), [2D Straight Skeleton and Polygon Offsetting](https://doc.cgal.org/6.0/Manual/packages.html#PkgStraightSkeleton2), [Surface Mesh](https://doc.cgal.org/6.0/Manual/packages.html#PkgSurfaceMesh), [2D Triangulations](https://doc.cgal.org/6.0/Manual/packages.html#PkgTriangulation2), [3D Triangulations](https://doc.cgal.org/6.0/Manual/packages.html#PkgTriangulation3), [2D Voronoi Diagrams](https://doc.cgal.org/6.0/Manual/packages.html#PkgVoronoiDiagram2), and more. The most simple use case of the basic viewer is the call of the global `CGAL::draw()` function. There is one such `draw()` function for each CGAL package that has a basic viewer. Such a call opens an interactive window showing the given model and allowing to navigate in the scene, show or hide some specific cells, show the interior of the model if any, etc. The `Basic_viewer` is based on Qt6.

### [Polygon Repair](https://doc.cgal.org/6.0/Manual/packages.html#PkgPolygonRepair) (new package)

- This package provides algorithms to repair 2D polygons, polygons with holes, and multipolygons with holes, by selecting faces of the arrangement of the input using the odd-even heuristic.

### [2D and 3D Linear Geometry Kernel](https://doc.cgal.org/6.0/Manual/packages.html#PkgKernel23)

- **Breaking change**: Replaced all instances of `boost::variant` with `std::variant` in the intersection functions.
- **Breaking change**: Replaced all instances of `boost::optional` with `std::optional` in the intersection functions.

### [3D Polyhedral Surface](https://doc.cgal.org/6.0/Manual/packages.html#PkgPolyhedron)

- The demo of this package, also known as “Polyhedron Demo” has been renamed “CGAL Lab” and moved to its own package (“Lab”).

### [2D and 3D Fast Intersection and Distance Computation (AABB Tree)](https://doc.cgal.org/6.0/Manual/packages.html#PkgAABBTree)

- The AABB tree can now be used with 2D or 3D primitives:
  - The concepts `AABBGeomTraits` and `AABBRayIntersectionGeomTraits` have been replaced by [`AABBGeomTraits_3`](https://doc.cgal.org/6.0/AABB_tree/classAABBGeomTraits__3.html) and by [`AABBRayIntersectionGeomTraits_3`](https://doc.cgal.org/6.0/AABB_tree/classAABBRayIntersectionGeomTraits__3.html), respectively.
  - The concepts [`AABBGeomTraits_2`](https://doc.cgal.org/6.0/AABB_tree/classAABBGeomTraits__2.html) and [`AABBRayIntersectionGeomTraits_2`](https://doc.cgal.org/6.0/AABB_tree/classAABBRayIntersectionGeomTraits__2.html) have been introduced, as the 2D counterparts.
  - The class [`CGAL::AABB_traits`](https://doc.cgal.org/6.0/AABB_tree/group__PkgAABBTreeRef.html#ga764f0fc59c96355877536810aa1aca5b) is deprecated and replaced by [`CGAL::AABB_traits_3`](https://doc.cgal.org/6.0/AABB_tree/classCGAL_1_1AABB__traits__3.html).
  - The class [`CGAL::AABB_traits_2`](https://doc.cgal.org/6.0/AABB_tree/classCGAL_1_1AABB__traits__2.html) is introduced as the 2D counterpart.
  - The class [`CGAL::AABB_segment_primitive`](https://doc.cgal.org/6.0/AABB_tree/group__PkgAABBTreeRef.html#gad0acfd5c4a3c081b7570cc6bd4594c8d) has been deprecated and replaced by the class [`CGAL::AABB_segment_primitive_3`](https://doc.cgal.org/6.0/AABB_tree/classCGAL_1_1AABB__segment__primitive__3.html).
  - The class [`CGAL::AABB_triangle_primitive`](https://doc.cgal.org/6.0/AABB_tree/group__PkgAABBTreeRef.html#ga54a56f01dc8024624f7d83ee0a01add0) has been deprecated and replaced by the class [`CGAL::AABB_triangle_primitive_3`](https://doc.cgal.org/6.0/AABB_tree/classCGAL_1_1AABB__triangle__primitive__3.html).
  - The following 2D primitive classes have been added: [`CGAL::AABB_segment_primitive_2`](https://doc.cgal.org/6.0/AABB_tree/classCGAL_1_1AABB__segment__primitive__2.html), [`CGAL::AABB_polyline_segment_primitive_2`](https://doc.cgal.org/6.0/AABB_tree/classCGAL_1_1AABB__polyline__segment__primitive__2.html), [`CGAL::AABB_triangle_primitive_2`](https://doc.cgal.org/6.0/AABB_tree/classCGAL_1_1AABB__triangle__primitive__2.html), [`CGAL::AABB_indexed_triangle_primitive_2`](https://doc.cgal.org/6.0/AABB_tree/classCGAL_1_1AABB__indexed__triangle__primitive__2.html).
- **Breaking change**: The concept [`AABBTraits`](https://doc.cgal.org/6.0/AABB_tree/classAABBTraits.html) now refines the concept [`SearchTraits`](https://doc.cgal.org/6.0/Spatial_searching/classSearchTraits.html).
- **Breaking change**: Replaced all instances of `boost::optional` with `std::optional`.

### [2D Arrangements](https://doc.cgal.org/6.0/Manual/packages.html#PkgArrangementOnSurface2)

- **Breaking change**: Replaced all instances of `boost::variant` with `std::variant`.
- **Breaking change**: The type of the result of point location queries has been changed to `std::variant`. Support for the old macro `CGAL_ARR_POINT_LOCATION_VERSION` has been removed.
- **Breaking change**: Eliminated the error-prone C-type casting that was used to define observers. In general, backward compatibility was maintained; however, the class template [`CGAL::Arr_observer`](https://doc.cgal.org/6.0/Arrangement_on_surface_2/group__PkgArrangementOnSurface2Ref.html#ga8019f986f5469920136c4b92290b7b1b) has been replaced by an alias template. (The class `CGAL::Arr_observer` was renamed to [`CGAL::Aos_observer`](https://doc.cgal.org/6.0/Arrangement_on_surface_2/classCGAL_1_1Aos__observer.html)).
- Introduced [`Arr_dcel`](https://doc.cgal.org/6.0/Arrangement_on_surface_2/classCGAL_1_1Arr__dcel.html), which essentially replaces the former `CGAL::Arr_default_dcel`. Backward compatibility was maintained by the introduction of the alias template [`CGAL::Arr_default_dcel`](https://doc.cgal.org/6.0/Arrangement_on_surface_2/group__PkgArrangementOnSurface2DCEL.html#gaf9635869a3794a46d7dcfce63d7de2a6). `CGAL::Arr_dcel`, as opposed to the former `CGAL::Arr_default_dcel` is templated (in addition to the geometry traits) by `Vertex`, `Halfedge`, and `Face` template parameters, and they have default type values. All this enables the layered extension of DCEL records.
- Fixed a bug in the zone construction code applied to arrangements of geodesic arcs on a sphere, when inserting an arc that lies on the identification curve.
- Introduced a new interactive program that demonstrates 2D arrangements embedded on the sphere called `earth`. The program (i) reads a database of all administrative boundaries of the countries in the world, (ii) displays the globe with all countries and land covered by water (which is land not covered by countries) on a window, and (ii) enables interaction with the user.

### [3D Envelopes](https://doc.cgal.org/6.0/Manual/packages.html#PkgEnvelope3)

- **Breaking change**: [`Construct_projected_boundary_2`](https://doc.cgal.org/6.0/Envelope_3/classEnvelopeTraits__3.html#ac7b8f72870f0572834a0a3de62c67bc1) in [`EnvelopeTraits_3`](https://doc.cgal.org/6.0/Envelope_3/classEnvelopeTraits__3.html) now uses `std::variant` instead of [`CGAL::Object`](https://doc.cgal.org/6.0/STL_Extension/classCGAL_1_1Object.html).
- Passed the base class of [`CGAL::Env_plane_traits_3`](https://doc.cgal.org/6.0/Envelope_3/classCGAL_1_1Env__plane__traits__3.html) as a template parameter with a default value (being the 2D arrangement linear traits). Similarly, passed the base class of `CGAL::Env_triangle_traits_3` as a template parameter with a default value (being the 2D arrangement segment traits).

### [Combinatorial Maps](https://doc.cgal.org/6.0/Manual/packages.html#PkgCombinatorialMaps) and [Generalized Maps](https://doc.cgal.org/6.0/Manual/packages.html#PkgGeneralizedMaps)

- Added the function [`insert_cell_1_between_two_cells_2()`](https://doc.cgal.org/6.0/Combinatorial_map/classGenericMap.html#aa29570a0812094c7876e24a228373f12) to the [`GenericMap`](https://doc.cgal.org/6.0/Combinatorial_map/classGenericMap.html) concept, which enables users to insert an edge between two different faces in order to create faces with holes.

- Added new meshing criterion `edge_distance`, an upper bound for the distance from the edge to the 1D feature.

- **Breaking change**: the concept `MeshEdgeCriteria_3` was modified to include the new meshing criterion `edge_distance`.

### [Quadtrees, Octrees, and Orthtrees](https://doc.cgal.org/6.0/Manual/packages.html#PkgOrthtree)

- **Breaking change**:
  - Node splitting behavior and per-node data are now customizable via the Traits class.
  - Nodes are now stored as a property map, with properties of each node accessed by index.
  - Nearest neighbors functions only work for Orthtrees which provide the necessary functionality.

### [CGAL and the Boost Graph Library (BGL)](https://doc.cgal.org/6.0/Manual/packages.html#PkgBGL)

- Added the function [`CGAL::remove_all_elements()`](https://doc.cgal.org/6.0/BGL/group__PkgBGLHelperFct.html#gac7e199820c95ed1fc6ab536750749358), which removes vertices, halfedges, and faces without collecting garbage and without removing properties.
- [Dynamic property maps](https://doc.cgal.org/6.0/BGL/group__PkgBGLPropertiesDynamic.html) can now have a default value.
- The class [`CGAL::Face_filtered_graph`](https://doc.cgal.org/6.0/BGL/structCGAL_1_1Face__filtered__graph.html) now supports patch IDs of any type and not just `faces_size_type`. The only requirement is that the type is hashable.

### [Polygon Mesh Processing](https://doc.cgal.org/6.0/Manual/packages.html#PkgPolygonMeshProcessing)

- Added the function [`CGAL::Polygon_mesh_processing::autorefine_triangle_soup()`](https://doc.cgal.org/6.0/Polygon_mesh_processing/group__PMP__corefinement__grp.html#gaec85370aa0b2acc0919e5f8406cfb74c), which can be used to refine a soup of triangles such that no pair of triangles intersects in their interiors. Also added, the function [`CGAL::Polygon_mesh_processing::autorefine()`](https://doc.cgal.org/6.0/Polygon_mesh_processing/group__PMP__corefinement__grp.html#ga3e3a0a82b6c04bdc3a6c070e8da4aed5) operating directly on a triangle mesh and updating it using the aforementioned function on a triangle soup.
- Added the class [`CGAL::Corefinement::Non_manifold_output_visitor`](https://doc.cgal.org/6.0/Polygon_mesh_processing/structCGAL_1_1Polygon__mesh__processing_1_1Corefinement_1_1Non__manifold__output__visitor.html), which can be used in corefinement based functions to deal with non-manifold outputs.
- Added the option to use a variable sizing field for [`CGAL::Polygon_mesh_processing::isotropic_remeshing()`](https://doc.cgal.org/6.0/Polygon_mesh_processing/group__PMP__meshing__grp.html#ga66cb01cf228ed22f0a2a474cfa2aeb3f), and a sizing function based on a measure of local curvature for adaptive remeshing.
- Added the function [`CGAL::Polygon_mesh_processing::interpolated_corrected_curvatures()`](https://doc.cgal.org/6.0/Polygon_mesh_processing/group__PMP__corrected__curvatures__grp.html#ga22665c9ce92aaedab07df1b05f20bdb2) which can be used to compute the mean and Gaussian curvatures, as well as the principal curvature and directions.
- Added the function [`CGAL::Polygon_mesh_processing::refine_mesh_at_isolevel()`](https://doc.cgal.org/6.0/Polygon_mesh_processing/group__PkgPolygonMeshProcessingRef.html#ga396505d5a60b5f6d29792b214fa59352) which can be used to refine a polygon mesh along an isocurve.
- Added the function [`CGAL::Polygon_mesh_processing::add_bbox()`](https://doc.cgal.org/6.0/Polygon_mesh_processing/group__PkgPolygonMeshProcessingRef.html#gabaf98d2fd9ae599ff1f3a5a6cde79cf3), which enables adding a tight or extended, triangulated or not, bounding box to a face graph.

### [2D Triangulations](https://doc.cgal.org/6.0/Manual/packages.html#PkgTriangulation2)

- **Breaking change**: the concept [`TriangulationTraits_2`](https://doc.cgal.org/6.0/Triangulation_2/classTriangulationTraits__2.html) now requires an additional functor `Compare_xy_2`.

### [3D Triangulations](https://doc.cgal.org/6.0/Manual/packages.html#PkgTriangulation3)

- Added three member functions [`vertices()`](https://doc.cgal.org/6.0/Triangulation_3/classCGAL_1_1Triangulation__3.html#a02faf334255e1ca8caa1a6f412533759) to the class [`CGAL::Triangulation_3`](https://doc.cgal.org/6.0/Triangulation_3/classCGAL_1_1Triangulation__3.html). Each of them returns an array containing the vertices of the given triangulation simplex.

### [dD Triangulations](https://doc.cgal.org/6.0/Manual/packages.html#PkgTriangulations)

- **Breaking change**: `CGAL::TDS_full_cell_mirror_storage_policy` is now unsupported in dimension larger than 127.
- **Breaking change**: Inserting multiple unweighted points in the same position now keeps the first one, instead of switching to the latest. This only affects custom point types where not all points in the same position are equivalent.

### [Tetrahedral Remeshing](https://doc.cgal.org/6.0/Manual/packages.html#PkgTetrahedralRemeshing)

- Added a sizing field as new parameter of [`CGAL::tetrahedral_isotropic_remeshing()`](https://doc.cgal.org/6.0/Tetrahedral_remeshing/group__PkgTetrahedralRemeshingRef.html#ga263775c52eeb483a86a16aeb9eb31af0), which can be used to perform non-uniform and adaptive remeshing.
- **Breaking change**: The template parameters of [`CGAL::Tetrahedral_remeshing::Remeshing_cell_base_3`](https://doc.cgal.org/6.0/Tetrahedral_remeshing/classCGAL_1_1Tetrahedral__remeshing_1_1Remeshing__cell__base__3.html) have been modified, reverting changes introduced in CGAL 5.6.
- **Breaking change**: The vertex base of [`CGAL::Tetrahedral_remeshing::Remeshing_vertex_base_3`](https://doc.cgal.org/6.0/Tetrahedral_remeshing/classCGAL_1_1Tetrahedral__remeshing_1_1Remeshing__vertex__base__3.html) must now be a model of the concept [`SimplicialMeshVertexBase_3`](https://doc.cgal.org/6.0/SMDS_3/classSimplicialMeshVertexBase__3.html) (and not only [`TriangulationVertexBase_3`](https://doc.cgal.org/6.0/Triangulation_3/classTriangulationVertexBase__3.html)).

### [3D Simplicial Mesh Data Structure](https://doc.cgal.org/6.0/Manual/packages.html#PkgSMDS3)

- **Breaking change**: The template parameters of [`CGAL::Simplicial_mesh_cell_base_3`](https://doc.cgal.org/6.0/SMDS_3/classCGAL_1_1Simplicial__mesh__cell__base__3.html) have been modified to enable passing a geometric traits and a custom cell base class.

### [3D Mesh Generation](https://doc.cgal.org/6.0/Manual/packages.html#PkgMesh3)

- **Breaking change**: Removed the concept `TriangleAccessor`, the template parameter `TriangleAccessor`, as well as the class `Triangle_accessor`. These were no longer used for several releases.
- **Breaking change**: Removed the class templates `CGAL::Gray_image_mesh_domain_3`, `CGAL::Implicit_mesh_domain_3`, and `CGAL::Labeled_image_mesh_domain_3`, which were deprecated since CGAL-4.13.

### [3D Surface Mesh Generation](https://doc.cgal.org/6.0/Manual/packages.html#PkgSurfaceMesher3)

- This package is deprecated and the package [3D Mesh Generation](https://doc.cgal.org/6.0/Manual/packages.html#PkgMesh3) should be used instead.

### [Surface Mesh Parameterization](https://doc.cgal.org/6.0/Manual/packages.html#PkgSurfaceMeshParameterization)

- **Breaking change**: The method [`CGAL::Surface_mesh_parameterization::LSCM_parameterizer_3`](https://doc.cgal.org/6.0/Surface_mesh_parameterization/classCGAL_1_1Surface__mesh__parameterization_1_1LSCM__parameterizer__3.html) now requires the Eigen library.
- **Breaking change**: CGAL no longer ships its own version of OpenNL.

### [Surface Mesh](https://doc.cgal.org/6.0/Manual/packages.html#PkgSurfaceMesh)

- **Breaking change**: The return type of [`CGAL::Surface_mesh::property_map()`](https://doc.cgal.org/6.0/Surface_mesh/classCGAL_1_1Surface__mesh.html#afc99c7ea179dc1c21a2ab59ed183184a) has been changed to `std::optional`.

### [3D Point Set](https://doc.cgal.org/6.0/Manual/packages.html#PkgPointSet3)

- **Breaking change**: The return type of [`CGAL::Point_set_3::property_map()`](https://doc.cgal.org/6.0/Point_set_3/classCGAL_1_1Point__set__3.html#a571ecc603cd32d78c7effaf86fe120ad) has been changed to `std::optional`.

### [Shape Detection](https://doc.cgal.org/6.0/Manual/packages.html#PkgShapeDetection)

- **Breaking change**: Replaced all instances of `boost::shared_ptr` with `std::shared_ptr`.

### [2D Straight Skeleton and Polygon Offsetting](https://doc.cgal.org/6.0/Manual/packages.html#PkgStraightSkeleton2)

- **Breaking change**: Replaced all instances of `boost::shared_ptr` with `std::shared_ptr`.
- **Breaking change**: Replaced all instances of `boost::optional` with `std::optional`.
