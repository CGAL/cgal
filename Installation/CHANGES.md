# Release History

## [Release 6.0.1](https://github.com/CGAL/cgal/releases/tag/v6.0.1)

### [Poisson Surface Reconstruction](https://doc.cgal.org/6.0.1/Manual/packages.html#PkgPoissonSurfaceReconstruction3)
-   Made the implicit function thread-safe so that the parallel version of `make_mesh_3()` can be used.

## [Release 6.0](https://github.com/CGAL/cgal/releases/tag/v6.0)

Release date: September 2024

### General Changes

- CGAL 6.0 is the first release of CGAL that requires a C++ compiler
  with the support of C++17 or later. The new list of supported compilers is:
  - Visual C++ 15.9, 16.10, 17.0 (from Visual Studio 2017, 2019 and 2022) or later
  - Gnu g++ 11.4.0 or later (on Linux or macOS)
  - LLVM Clang version 15.0.7 or later (on Linux)
  - Apple Clang compiler versions 10.0.1, 12.0.5, and 15.0.0 (on macOS)
- The minimal supported version of Boost is now 1.72.0.
- GMP/MPFR are no longer mandatory to use CGAL, [Boost.Multiprecision](https://www.boost.org/doc/libs/1_72_0/libs/multiprecision/doc/html/index.html).
  can be used instead.
- The CGAL `Core` library is no longer based on GMP, but on
  [Boost.Multiprecision](https://www.boost.org/doc/libs/1_72_0/libs/multiprecision/doc/html/index.html).
  Either GMP backend or Boost backend can be used.
- All demos are now based on Qt6.
- **Breaking change**: The CMake file `UseCGAL.cmake` has been removed from CGAL.
  Usages of the CMake variables `${CGAL_USE_FILE}` and `${CGAL_LIBRARIES}` must be replaced
  by a link to the imported target `CGAL::CGAL`, for example:
  `target_link_library(your_target PRIVATE CGAL::CGAL)`.

### [Kinetic Space Partition](https://doc.cgal.org/6.0/Manual/packages.html#PkgKineticSpacePartition) (new package)

-   This package implements kinetic space partition: based on a set of planar input shapes,
    the bounding box of the input data is split into convex volumes. The complexity of the partition
    can be adjusted with a single parameter.

### [Kinetic Surface Reconstruction](https://doc.cgal.org/6.0/Manual/packages.html#PkgKineticSurfaceReconstruction) (new package)

-   The package implements a piece-wise planar surface reconstruction pipeline from point clouds
    combining methods from the [Shape Detection](https://doc.cgal.org/6.0/Manual/packages.html#PkgShapeDetection),
    [Shape Regularization](https://doc.cgal.org/6.0/Manual/packages.html#PkgShapeRegularization)
    and [Kinetic Shape Partition](https://doc.cgal.org/6.0/Manual/packages.html#PkgKineticSpacePartition) packages
    and graph-cut to reconstruct surfaces from point clouds.

### [Basic Viewer](https://doc.cgal.org/6.0/Basic_viewer/index.html#Chapter_Basic_viewer) (new package)

-   The basic viewer package provides interactive visualization for most CGAL packages,
    such as [2D Arrangements](https://doc.cgal.org/6.0/Manual/packages.html#PkgArrangementOnSurface2),
    [2D Regularized Boolean Set-Operations](https://doc.cgal.org/6.0/Manual/packages.html#PkgBooleanSetOperations2),
    [Linear Cell Complex](https://doc.cgal.org/6.0/Manual/packages.html#PkgLinearCellComplex),
    [3D Boolean Operations on Nef Polyhedra](https://doc.cgal.org/6.0/Manual/packages.html#PkgNef3),
    [2D Periodic Triangulations](https://doc.cgal.org/6.0/Manual/packages.html#PkgPeriodic2Triangulation2),
    [3D Point Set](https://doc.cgal.org/6.0/Manual/packages.html#PkgPointSet3),
    [2D Polygons](https://doc.cgal.org/6.0/Manual/packages.html#PkgPolygon2),
    [3D Polyhedral Surface](https://doc.cgal.org/6.0/Manual/packages.html#PkgPolyhedron),
    [2D Straight Skeleton and Polygon Offsetting](https://doc.cgal.org/6.0/Manual/packages.html#PkgStraightSkeleton2),
    [Surface Mesh](https://doc.cgal.org/6.0/Manual/packages.html#PkgSurfaceMesh),
    [2D Triangulations](https://doc.cgal.org/6.0/Manual/packages.html#PkgTriangulation2),
    [3D Triangulations](https://doc.cgal.org/6.0/Manual/packages.html#PkgTriangulation3),
    [2D Voronoi Diagrams](https://doc.cgal.org/6.0/Manual/packages.html#PkgVoronoiDiagram2),
    and more.
    The most simple use case of the basic viewer is the call of the global `CGAL::draw()` function.
    There is one such `draw()` function for each CGAL package that has a basic viewer. Such a call opens
    an interactive window showing the given model and allowing to navigate in the scene,
    show or hide some specific cells, show the interior of the model if any, etc.
    The `Basic_viewer` is based on Qt6.

### [Polygon Repair](https://doc.cgal.org/6.0/Manual/packages.html#PkgPolygonRepair) (new package)

-   This package provides algorithms to repair 2D polygons, polygons with holes,
    and multipolygons with holes, by selecting faces of the arrangement of the input
    using the odd-even heuristic.

### [2D and 3D Linear Geometry Kernel](https://doc.cgal.org/6.0/Manual/packages.html#PkgKernel23)

-   **Breaking change**: Replaced all instances of `boost::variant` with `std::variant`
    in the intersection functions.
-   **Breaking change**: Replaced all instances of `boost::optional` with `std::optional`
    in the intersection functions.

### [3D Polyhedral Surface](https://doc.cgal.org/6.0/Manual/packages.html#PkgPolyhedron)

-   The demo of this package, also known as "Polyhedron Demo" has been renamed "CGAL Lab"
    and moved to its own package ("Lab").

### [2D and 3D Fast Intersection and Distance Computation (AABB Tree)](https://doc.cgal.org/6.0/Manual/packages.html#PkgAABBTree)

- The AABB tree can now be used with 2D or 3D primitives:
  - The concepts `AABBGeomTraits` and `AABBRayIntersectionGeomTraits`
    have been replaced by [`AABBGeomTraits_3`](https://doc.cgal.org/6.0/AABB_tree/classAABBGeomTraits__3.html)
    and by [`AABBRayIntersectionGeomTraits_3`](https://doc.cgal.org/6.0/AABB_tree/classAABBRayIntersectionGeomTraits__3.html),
    respectively.
  - The concepts [`AABBGeomTraits_2`](https://doc.cgal.org/6.0/AABB_tree/classAABBGeomTraits__2.html)
    and [`AABBRayIntersectionGeomTraits_2`](https://doc.cgal.org/6.0/AABB_tree/classAABBRayIntersectionGeomTraits__2.html)
    have been introduced, as the 2D counterparts.
  - The class [`CGAL::AABB_traits`](https://doc.cgal.org/6.0/AABB_tree/group__PkgAABBTreeRef.html#ga764f0fc59c96355877536810aa1aca5b)
    is deprecated and replaced by [`CGAL::AABB_traits_3`](https://doc.cgal.org/6.0/AABB_tree/classCGAL_1_1AABB__traits__3.html).
  - The class [`CGAL::AABB_traits_2`](https://doc.cgal.org/6.0/AABB_tree/classCGAL_1_1AABB__traits__2.html) is introduced as the 2D counterpart.
  - The class [`CGAL::AABB_segment_primitive`](https://doc.cgal.org/6.0/AABB_tree/group__PkgAABBTreeRef.html#gad0acfd5c4a3c081b7570cc6bd4594c8d)
    has been deprecated and replaced by the class [`CGAL::AABB_segment_primitive_3`](https://doc.cgal.org/6.0/AABB_tree/classCGAL_1_1AABB__segment__primitive__3.html).
  - The class [`CGAL::AABB_triangle_primitive`](https://doc.cgal.org/6.0/AABB_tree/group__PkgAABBTreeRef.html#ga54a56f01dc8024624f7d83ee0a01add0)
    has been deprecated and replaced by the class [`CGAL::AABB_triangle_primitive_3`](https://doc.cgal.org/6.0/AABB_tree/classCGAL_1_1AABB__triangle__primitive__3.html).
  - The following 2D primitive classes have been added:
    [`CGAL::AABB_segment_primitive_2`](https://doc.cgal.org/6.0/AABB_tree/classCGAL_1_1AABB__segment__primitive__2.html),
    [`CGAL::AABB_polyline_segment_primitive_2`](https://doc.cgal.org/6.0/AABB_tree/classCGAL_1_1AABB__polyline__segment__primitive__2.html),
    [`CGAL::AABB_triangle_primitive_2`](https://doc.cgal.org/6.0/AABB_tree/classCGAL_1_1AABB__triangle__primitive__2.html),
    [`CGAL::AABB_indexed_triangle_primitive_2`](https://doc.cgal.org/6.0/AABB_tree/classCGAL_1_1AABB__indexed__triangle__primitive__2.html).
- **Breaking change**: The concept [`AABBTraits`](https://doc.cgal.org/6.0/AABB_tree/classAABBTraits.html)
    now refines the concept [`SearchTraits`](https://doc.cgal.org/6.0/Spatial_searching/classSearchTraits.html).
- **Breaking change**: Replaced all instances of `boost::optional` with `std::optional`.

### [2D Arrangements](https://doc.cgal.org/6.0/Manual/packages.html#PkgArrangementOnSurface2)

-   **Breaking change**: Replaced all instances of `boost::variant` with `std::variant`.
-   **Breaking change**: The type of the result of point location queries has been changed to
    `std::variant`. Support for the old macro `CGAL_ARR_POINT_LOCATION_VERSION`
    has been removed.
-   **Breaking change**: Eliminated the error-prone C-type casting that was used to define observers.
    In general, backward compatibility was maintained; however, the class template
    [`CGAL::Arr_observer`](https://doc.cgal.org/6.0/Arrangement_on_surface_2/group__PkgArrangementOnSurface2Ref.html#ga8019f986f5469920136c4b92290b7b1b)
    has been replaced by an alias template. (The class `CGAL::Arr_observer`
    was renamed to [`CGAL::Aos_observer`](https://doc.cgal.org/6.0/Arrangement_on_surface_2/classCGAL_1_1Aos__observer.html)).
-   Introduced [`Arr_dcel`](https://doc.cgal.org/6.0/Arrangement_on_surface_2/classCGAL_1_1Arr__dcel.html),
    which essentially replaces the former `CGAL::Arr_default_dcel`.
    Backward compatibility was maintained by the introduction of the alias template
    [`CGAL::Arr_default_dcel`](https://doc.cgal.org/6.0/Arrangement_on_surface_2/group__PkgArrangementOnSurface2DCEL.html#gaf9635869a3794a46d7dcfce63d7de2a6).
    `CGAL::Arr_dcel`, as opposed to the former `CGAL::Arr_default_dcel` is templated
    (in addition to the geometry traits) by `Vertex`, `Halfedge`, and `Face` template parameters,
    and they have default type values. All this enables the layered extension of DCEL records.
-   Fixed a bug in the zone construction code applied to arrangements of geodesic arcs on a sphere,
    when inserting an arc that lies on the identification curve.
-   Introduced a new interactive program that demonstrates 2D arrangements embedded on the sphere
    called `earth`. The program (i) reads a database of all administrative boundaries of the countries
    in the world, (ii) displays the globe with all countries and land covered by water (which is land
    not covered by countries) on a window, and (ii) enables interaction with the user.

### [3D Envelopes](https://doc.cgal.org/6.0/Manual/packages.html#PkgEnvelope3)

-   **Breaking change**: [`Construct_projected_boundary_2`](https://doc.cgal.org/6.0/Envelope_3/classEnvelopeTraits__3.html#ac7b8f72870f0572834a0a3de62c67bc1)
    in [`EnvelopeTraits_3`](https://doc.cgal.org/6.0/Envelope_3/classEnvelopeTraits__3.html)
    now uses `std::variant` instead of [`CGAL::Object`](https://doc.cgal.org/6.0/STL_Extension/classCGAL_1_1Object.html).
-    Passed the base class of [`CGAL::Env_plane_traits_3`](https://doc.cgal.org/6.0/Envelope_3/classCGAL_1_1Env__plane__traits__3.html)
    as a template parameter with a default value (being the 2D arrangement linear traits).
    Similarly, passed the base class of `CGAL::Env_triangle_traits_3` as a template parameter
    with a default value (being the 2D arrangement segment traits).

### [Combinatorial Maps](https://doc.cgal.org/6.0/Manual/packages.html#PkgCombinatorialMaps) and [Generalized Maps](https://doc.cgal.org/6.0/Manual/packages.html#PkgGeneralizedMaps)

-    Added the function [`insert_cell_1_between_two_cells_2()`](https://doc.cgal.org/6.0/Combinatorial_map/classGenericMap.html#aa29570a0812094c7876e24a228373f12)
    to the [`GenericMap`](https://doc.cgal.org/6.0/Combinatorial_map/classGenericMap.html)
    concept, which enables users to insert an edge between two different faces in order to create faces with holes.

### [Quadtrees, Octrees, and Orthtrees](https://doc.cgal.org/6.0/Manual/packages.html#PkgOrthtree)

- **Breaking change**:
  - Node splitting behavior and per-node data are now customizable via the Traits class.
  - Nodes are now stored as a property map, with properties of each node accessed by index.
  - Nearest neighbors functions only work for Orthtrees which provide the necessary functionality.

### [CGAL and the Boost Graph Library (BGL)](https://doc.cgal.org/6.0/Manual/packages.html#PkgBGL)

-   Added the function [`CGAL::remove_all_elements()`](https://doc.cgal.org/6.0/BGL/group__PkgBGLHelperFct.html#gac7e199820c95ed1fc6ab536750749358),
    which removes vertices, halfedges, and faces without collecting garbage and without removing properties.
-   [Dynamic property maps](https://doc.cgal.org/6.0/BGL/group__PkgBGLPropertiesDynamic.html)
    can now have a default value.
-   The class [`CGAL::Face_filtered_graph`](https://doc.cgal.org/6.0/BGL/structCGAL_1_1Face__filtered__graph.html)
    now supports patch IDs of any type and not just `faces_size_type`. The only requirement is that
    the type is hashable.

### [Polygon Mesh Processing](https://doc.cgal.org/6.0/Manual/packages.html#PkgPolygonMeshProcessing)

-   Added the function [`CGAL::Polygon_mesh_processing::autorefine_triangle_soup()`](https://doc.cgal.org/6.0/Polygon_mesh_processing/group__PMP__corefinement__grp.html#gaec85370aa0b2acc0919e5f8406cfb74c),
    which can be used to refine a soup of triangles such that no pair of triangles intersects
    in their interiors. Also added, the function [`CGAL::Polygon_mesh_processing::autorefine()`](https://doc.cgal.org/6.0/Polygon_mesh_processing/group__PMP__corefinement__grp.html#ga3e3a0a82b6c04bdc3a6c070e8da4aed5)
    operating directly on a triangle mesh and updating it using the aforementioned function on a triangle soup.
-   Added the class [`CGAL::Corefinement::Non_manifold_output_visitor`](https://doc.cgal.org/6.0/Polygon_mesh_processing/structCGAL_1_1Polygon__mesh__processing_1_1Corefinement_1_1Non__manifold__output__visitor.html),
    which can be used in corefinement based functions to deal with non-manifold outputs.
-   Added the option to use a variable sizing field for [`CGAL::Polygon_mesh_processing::isotropic_remeshing()`](https://doc.cgal.org/6.0/Polygon_mesh_processing/group__PMP__meshing__grp.html#ga66cb01cf228ed22f0a2a474cfa2aeb3f),
    and a sizing function based on a measure of local curvature for adaptive remeshing.
-   Added the function [`CGAL::Polygon_mesh_processing::interpolated_corrected_curvatures()`](https://doc.cgal.org/6.0/Polygon_mesh_processing/group__PMP__corrected__curvatures__grp.html#ga22665c9ce92aaedab07df1b05f20bdb2)
    which can be used to compute the mean and Gaussian curvatures, as well as the principal curvature
    and directions.
-   Added the function [`CGAL::Polygon_mesh_processing::refine_mesh_at_isolevel()`](https://doc.cgal.org/6.0/Polygon_mesh_processing/group__PkgPolygonMeshProcessingRef.html#ga396505d5a60b5f6d29792b214fa59352)
    which can be used to refine a polygon mesh along an isocurve.
-   Added the function [`CGAL::Polygon_mesh_processing::add_bbox()`](https://doc.cgal.org/6.0/Polygon_mesh_processing/group__PkgPolygonMeshProcessingRef.html#gabaf98d2fd9ae599ff1f3a5a6cde79cf3),
    which enables adding a tight or extended, triangulated or not, bounding box to a face graph.

### [2D Triangulations](https://doc.cgal.org/6.0/Manual/packages.html#PkgTriangulation2)
-   **Breaking change**: the concept [`TriangulationTraits_2`](https://doc.cgal.org/6.0/Triangulation_2/classTriangulationTraits__2.html) now requires an additional functor `Compare_xy_2`.

### [3D Triangulations](https://doc.cgal.org/6.0/Manual/packages.html#PkgTriangulation3)

-   Added three member functions [`vertices()`](https://doc.cgal.org/6.0/Triangulation_3/classCGAL_1_1Triangulation__3.html#a02faf334255e1ca8caa1a6f412533759)
    to the class [`CGAL::Triangulation_3`](https://doc.cgal.org/6.0/Triangulation_3/classCGAL_1_1Triangulation__3.html).
    Each of them returns an array containing the vertices of the given triangulation simplex.

### [dD Triangulations](https://doc.cgal.org/6.0/Manual/packages.html#PkgTriangulations)

-   **Breaking change**: `CGAL::TDS_full_cell_mirror_storage_policy` is now unsupported in dimension larger than 127.
-   **Breaking change**: Inserting multiple unweighted points in the same
    position now keeps the first one, instead of switching to the latest. This
    only affects custom point types where not all points in the same position
    are equivalent.

### [Tetrahedral Remeshing](https://doc.cgal.org/6.0/Manual/packages.html#PkgTetrahedralRemeshing)

-   Added a sizing field as new parameter of [`CGAL::tetrahedral_isotropic_remeshing()`](https://doc.cgal.org/6.0/Tetrahedral_remeshing/group__PkgTetrahedralRemeshingRef.html#ga263775c52eeb483a86a16aeb9eb31af0),
    which can be used to perform non-uniform and adaptive remeshing.
-   **Breaking change**: The template parameters of
    [`CGAL::Tetrahedral_remeshing::Remeshing_cell_base_3`](https://doc.cgal.org/6.0/Tetrahedral_remeshing/classCGAL_1_1Tetrahedral__remeshing_1_1Remeshing__cell__base__3.html)
    have been modified, reverting changes introduced in CGAL 5.6.
-   **Breaking change**: The vertex base of
    [`CGAL::Tetrahedral_remeshing::Remeshing_vertex_base_3`](https://doc.cgal.org/6.0/Tetrahedral_remeshing/classCGAL_1_1Tetrahedral__remeshing_1_1Remeshing__vertex__base__3.html)
    must now be a model of the concept
    [`SimplicialMeshVertexBase_3`](https://doc.cgal.org/6.0/SMDS_3/classSimplicialMeshVertexBase__3.html)
    (and not only [`TriangulationVertexBase_3`](https://doc.cgal.org/6.0/Triangulation_3/classTriangulationVertexBase__3.html)).

### [3D Simplicial Mesh Data Structure](https://doc.cgal.org/6.0/Manual/packages.html#PkgSMDS3)

-   **Breaking change**: The template parameters of
    [`CGAL::Simplicial_mesh_cell_base_3`](https://doc.cgal.org/6.0/SMDS_3/classCGAL_1_1Simplicial__mesh__cell__base__3.html)
    have been modified to enable passing a geometric traits and a custom cell base class.

### [3D Mesh Generation](https://doc.cgal.org/6.0/Manual/packages.html#PkgMesh3)

-   **Breaking change**: Removed the concept `TriangleAccessor`, the template parameter `TriangleAccessor`,
    as well as the class `Triangle_accessor`. These were no longer used for several releases.
-   **Breaking change**: Removed the class templates `CGAL::Gray_image_mesh_domain_3`, `CGAL::Implicit_mesh_domain_3`,
    and `CGAL::Labeled_image_mesh_domain_3`, which were deprecated since CGAL-4.13.
-   Added new meshing criterion `edge_distance`, an upper bound for the distance from the edge to the 1D feature.
- **Breaking change**: the concept `MeshEdgeCriteria_3` was modified to include the new meshing criterion `edge_distance`.


### [3D Surface Mesh Generation](https://doc.cgal.org/6.0/Manual/packages.html#PkgSurfaceMesher3)

-   This package is deprecated and the package [3D Mesh Generation](https://doc.cgal.org/6.0/Manual/packages.html#PkgMesh3) should
    be used instead.

### [Surface Mesh Parameterization](https://doc.cgal.org/6.0/Manual/packages.html#PkgSurfaceMeshParameterization)

-   **Breaking change**: The method [`CGAL::Surface_mesh_parameterization::LSCM_parameterizer_3`](https://doc.cgal.org/6.0/Surface_mesh_parameterization/classCGAL_1_1Surface__mesh__parameterization_1_1LSCM__parameterizer__3.html)
    now requires the Eigen library.
-   **Breaking change**: CGAL no longer ships its own version of OpenNL.

### [Surface Mesh](https://doc.cgal.org/6.0/Manual/packages.html#PkgSurfaceMesh)

-   **Breaking change**: The return type of [`CGAL::Surface_mesh::property_map()`](https://doc.cgal.org/6.0/Surface_mesh/classCGAL_1_1Surface__mesh.html#afc99c7ea179dc1c21a2ab59ed183184a)
    has been changed to `std::optional`.

### [3D Point Set](https://doc.cgal.org/6.0/Manual/packages.html#PkgPointSet3)

-   **Breaking change**: The return type of [`CGAL::Point_set_3::property_map()`](https://doc.cgal.org/6.0/Point_set_3/classCGAL_1_1Point__set__3.html#a571ecc603cd32d78c7effaf86fe120ad)
    has been changed to `std::optional`.

### [Shape Detection](https://doc.cgal.org/6.0/Manual/packages.html#PkgShapeDetection)

-   **Breaking change**: Replaced all instances of `boost::shared_ptr` with `std::shared_ptr`.

### [2D Straight Skeleton and Polygon Offsetting](https://doc.cgal.org/6.0/Manual/packages.html#PkgStraightSkeleton2)

-   **Breaking change**: Replaced all instances of `boost::shared_ptr` with `std::shared_ptr`.
-   **Breaking change**: Replaced all instances of `boost::optional` with `std::optional`.


[Release 5.6](https://github.com/CGAL/cgal/releases/tag/v5.6)
-----------

Release date: July 2023

### General Changes

-   **Breaking change**: Package-specific assertions, preconditions, and postconditions (such as
    `CGAL_triangulation_assertion`) have been removed. Corresponding CGAL-wide versions (such as
    `CGAL_assertion`) should be used instead.

### [Shape Detection](https://doc.cgal.org/5.6/Manual/packages.html#PkgShapeDetection) (major changes)

-   **Breaking change**: The region growing part of the package have been reworked to fix design
    issues introduced with the handling of `FaceGraph` models. In particular, the notion of `Item`
    has been introduced to reference an element in the input range of elements. Region maps now
    operates on `Item` and no longer on the value type of the input range.
-   **Breaking change**: The method `update()` in the concept `RegionType` now returns a `Boolean`
    instead of `void`, that is used inside the class `Region_growing` for detecting if the input
    conditions for the new region are satisfied. This change affects only user-defined types of regions.
-   **Breaking change**: The constructors of all models used together with the region growing algorithm
    now enable users to provide parameters through the [named parameters](https://doc.cgal.org/5.6/BGL/group__bgl__namedparameters.html) mechanism.
-   All fitting classes in the region growing framework are now using better versions of the region
    conditions, more precise and faster, including the correct normal orientations.
-   Added new models of the concept `RegionType` for getting linear regions in a set of 2D and 3D
    segments and on 2D and 3D polylines.
-   Added the class `Polyline_graph` for extracting a set of polylines from a face graph, which splits
    this graph into a set of user-defined regions.
-   Added new shapes to the Region Growing algorithm on a point set: circles in 2D, spheres in 3D,
    and cylinders in 3D.

### [2D Straight Skeleton and Polygon Offsetting](https://doc.cgal.org/5.6/Manual/packages.html#PkgStraightSkeleton2) (major changes)
-   Added weighted straight skeletons: weighted straight skeletons are a generalization of
    straight skeletons. Contour edges are assigned a positive weight, which can be understood
    as assigning a speed to the wavefront spawned from the contour edge.
-   Added straight skeleton extrusion: this CGAL package now implements the extrusion of weighted
    straight skeletons of polygons with holes. The output is a closed, combinatorially 2-manifold
    surface triangle mesh.

    See also the [news entry](https://www.cgal.org/2023/05/09/improved_straight_skeleton/).

### [2D and 3D Linear Geometry Kernel](https://doc.cgal.org/5.6/Manual/packages.html#PkgKernel23)

-   Added the functor
    [`CompareAngle_3`](https://doc.cgal.org/5.6/Kernel_23/classKernel_1_1CompareAngle__3.html)
    to the concept
    [`Kernel`](https://doc.cgal.org/5.6/Kernel_23/classKernel.html) to compare
    an angle defined by three points to the cosinus of another angle.

### [Combinatorial Maps](https://doc.cgal.org/5.6/Manual/packages.html#PkgCombinatorialMaps), [Generalized Maps](https://doc.cgal.org/5.6/Manual/packages.html#PkgGeneralizedMaps), and [Linear Cell Complex](https://doc.cgal.org/5.6/Manual/packages.html#PkgLinearCellComplex)

-   Added a version that uses indices instead of handles as dart and attribute descriptors.
    As the indices are integers convertible from and to `std::size_t`, they can be used as index
    into vectors which store properties. To use the index version, `Use_index` must be defined
    and be equal to `CGAL::Tag_true` in the item class.

### [Linear Cell Complex](https://doc.cgal.org/5.6/Manual/packages.html#PkgLinearCellComplex)

-   Added the class
    [`Linear_cell_complex_incremental_builder_3`](https://doc.cgal.org/5.6/Linear_cell_complex/classCGAL_1_1Linear__cell__complex__incremental__builder__3.html).

### [2D Arrangements](https://doc.cgal.org/5.6/Manual/packages.html#PkgArrangementOnSurface2)

-   Introduced an overload function template, namely `draw(arr)`, that renders arrangements based
    on the `Basic_viewer_qt` class template. As of now, only 2D arrangements on the plane induced
    by (i) segments, (ii) conics, and (iii) circular arcs or (linear) segments are supported.
-   Improved the traits class template that handles conics, namely
    [`Arr_conic_traits_2`](https://doc.cgal.org/5.6/Arrangement_on_surface_2/classCGAL_1_1Arr__conic__traits__2.html).
    This includes the following:
    1. Fixed a couple of bugs and slightly optimized some functions.
    2. Introduced functionality that approximates conics with polylines. (This is used to draw conic curves.)
    3. **Breaking change**: Changed the interface to generate conic curves. In the past, curves where
    generated directly using the constructors of the conic and x-monotone conic constructs. Now,
    they are constructed via function objects provided by the traits. This eliminates the constructions
    of temporary kernels. The old functionality is obsolete, but still supported for a limited number
    of versions. It depends on a static member function of the traits. In a future version this function
    will no longer be static, implying that the old functionality will no longer be supported.
- Introduced functionality that approximates circular segments with polylines. (This is used to draw conic curves.)

### [Polygon Mesh Processing](https://doc.cgal.org/5.6/Manual/packages.html#PkgPolygonMeshProcessing)

-   Added functions
    [`CGAL::Polygon_mesh_processing::region_growing_of_planes_on_faces()`](https://doc.cgal.org/5.6/Polygon_mesh_processing/group__PkgPolygonMeshProcessingRef.html#ga50dcd2f6295f584d2e378b57290ae2af)
    and
    [`CGAL::Polygon_mesh_processing::detect_corners_of_regions()`](https://doc.cgal.org/5.6/Polygon_mesh_processing/group__PkgPolygonMeshProcessingRef.html#gac8e445730d718a2fc49604e865017d2e),
    which enable partitioning a mesh into planar regions using the region growing algorithm
    from the [Shape Detection](https://doc.cgal.org/5.6/Manual/packages.html#PkgShapeDetection) package.

-   Added the functions
    [`CGAL::Polygon_mesh_processing::remesh_planar_patches()`](https://doc.cgal.org/5.6/Polygon_mesh_processing/group__PMP__meshing__grp.html#ga7fca6fa2db94560ab6d32e6a77fc35b6)
    and
    [`CGAL::Polygon_mesh_processing::remesh_almost_planar_patches()`](https://doc.cgal.org/5.6/Polygon_mesh_processing/group__PMP__meshing__grp.html#ga0e6da479548199a5d82c3cf0ed36e8a0),
    which can be used to remesh patches of coplanar faces in a mesh.

-   Added the function
    [`CGAL::Polygon_mesh_processing::surface_Delaunay_remeshing()`](https://doc.cgal.org/5.6/Polygon_mesh_processing/group__PMP__meshing__grp.html#gaff62f9415d2fe96d1d3095351f156ced),
    which can be used to remesh a surface triangle mesh using the Delaunay refinement algorithm
    from the [3D Mesh Generation](https://doc.cgal.org/5.6/Manual/packages.html#PkgMesh3) package.

-   Added the function
    [`CGAL::Polygon_mesh_processing::remove_almost_degenerate_faces()`](https://doc.cgal.org/5.6/Polygon_mesh_processing/group__PMP__geometric__repair__grp.html#ga48008d2b66de8a68a7068f29db15dad6),
    which can be used to remove badly shaped triangles faces in a mesh.

-   Added the functions
    [`CGAL::Polygon_mesh_processing::does_triangle_soup_self_intersect()`](https://doc.cgal.org/5.6/Polygon_mesh_processing/group__PMP__intersection__grp.html#ga4909920dc48b8285e69feb845feb1e53)
    and
    [`CGAL::Polygon_mesh_processing::triangle_soup_self_intersections()`](https://doc.cgal.org/5.6/Polygon_mesh_processing/group__PMP__intersection__grp.html#ga1c5fee17bd0d92d5a2fba77ed94d4b4d)
    to identify and report self-intersections in a triangle soup, similarly to existing functions on triangle meshes.

-   Added the function
    [`CGAL::Polygon_mesh_processing::triangulate_polygons()`](https://doc.cgal.org/5.6/Polygon_mesh_processing/group__PMP__meshing__grp.html#ga8b7db6aa8c3e79526b594739ba926d82),
    which allows users to triangulate polygon soups.

-   Added a named parameter to
    [`CGAL::Polygon_mesh_processing::smooth_shape()`](https://doc.cgal.org/5.6/Polygon_mesh_processing/group__PMP__meshing__grp.html#ga57fa999abe8dc557003482444df2a189)
    to disable the scaling, which otherwise aims to compensate volume loss during smoothing.

-   Deprecated the overloads of functions
    [`CGAL::Polygon_mesh_processing::triangulate_hole()`](https://doc.cgal.org/5.6/Polygon_mesh_processing/group__PMP__hole__filling__grp.html#ga3abdf2d0558822e85f060966b69cae98),
    [`CGAL::Polygon_mesh_processing::triangulate_and_refine_hole()`](https://doc.cgal.org/5.6/Polygon_mesh_processing/group__PMP__hole__filling__grp.html#ga9868fac4d9dca77462ad7828bc99d8a1),
    and
    [`CGAL::Polygon_mesh_processing::triangulate_refine_and_fair_hole()`](https://doc.cgal.org/5.6/Polygon_mesh_processing/group__PMP__hole__filling__grp.html#ga18eac756a8f8e5d5f73e645fd4e26cad)
    which have output iterators for vertices and faces as parameter. They are replaced by overloads
    with two additional named parameters.

### [2D Convex Hulls](https://doc.cgal.org/5.6/Manual/packages.html#PkgConvexHull2)

-   **Breaking change**: The concept
    [`ConvexHullTraits_2`](https://doc.cgal.org/5.6/Convex_hull_2/classConvexHullTraits__2.html)
    no longer requires the functor `Less_signed_distance_to_line_2`, but requires the functor
    `Compare_signed_distance_to_line_2`
    instead.
-   The long-deprecated classes `Convex_hull_projective_xy_traits_2`, `Convex_hull_projective_xz_traits_2`,
    and `Convex_hull_projective_yz_traits_2` have been removed. Users should use
    [`Projection_traits_xy_3`](https://doc.cgal.org/5.6/Kernel_23/classCGAL_1_1Projection__traits__xy__3.html),
    [`Projection_traits_xz_3`](https://doc.cgal.org/5.6/Kernel_23/classCGAL_1_1Projection__traits__xz__3.html),
    and
    [`Projection_traits_yz_3`](https://doc.cgal.org/5.6/Kernel_23/classCGAL_1_1Projection__traits__yz__3.html)
    instead.

### [2D Triangulations](https://doc.cgal.org/5.6/Manual/packages.html#PkgTriangulation2)

-   Added the function
    [`CGAL::mark_domain_in_triangulation()`](https://doc.cgal.org/5.6/Triangulation_2/group__PkgTriangulation2Miscellaneous.html#ga0409755d0eb89100810230443a85e7eb)
    to mark faces connected with non-constrained edges as inside of the domain based on the nesting level.

### [2D Conforming Triangulations and Meshes](https://doc.cgal.org/5.6/Manual/packages.html#PkgMesh2)

-   Added new overloads to the function
    [`write_VTU()`](https://doc.cgal.org/5.6/Mesh_2/group__PkgMesh2IO.html),
    with property maps for specifying the domain.
-   Deprecated usage of boost parameters in favor of function named parameters in
    [`CGAL::lloyd_optimize_mesh_2()`](https://doc.cgal.org/5.6/Mesh_2/group__PkgMesh2Functions.html#gafeaf59d3fa014da287f8514913b38d05).
-   Deprecated two overloads of the function
    [`refine_Delaunay_mesh()`](https://doc.cgal.org/5.6/Mesh_2/group__PkgMesh2Functions.html),
    and replaced them with versions using function named parameters.

### [2D Hyperbolic Triangulations](https://doc.cgal.org/5.6/Manual/packages.html#PkgHyperbolicTriangulation2)

-   **Breaking change**: the concept
    [`HyperbolicTriangulationFaceBase_2`](https://doc.cgal.org/5.6/Hyperbolic_triangulation_2/classHyperbolicTriangulationFaceBase__2.html)
    has been modified to
    better reflect the triangulation's requirements and avoid a conflict with the requirements
    described by the concept `TriangulationDataStructure_2::Face`. The model
    [`CGAL::Hyperbolic_triangulation_face_base_2`](https://doc.cgal.org/5.6/Hyperbolic_triangulation_2/classCGAL_1_1Hyperbolic__triangulation__face__base__2.html)
    has been adapted correspondingly.

### [3D Simplicial Mesh Data Structure](https://doc.cgal.org/5.6/Manual/packages.html#PkgSMDS3) (new package)

-   This new package wraps all the existing code that deals with a
    [`MeshComplex_3InTriangulation_3`](https://doc.cgal.org/5.6/SMDS_3/classMeshComplex__3InTriangulation__3.html)
    to describe 3D simplicial meshes, and makes the data structure independent
    from the [tetrahedral mesh generation](https://doc.cgal.org/5.6/Manual/packages.html#PkgMesh3) package.

### [Tetrahedral Remeshing](https://doc.cgal.org/5.6/Manual/packages.html#PkgTetrahedralRemeshing)
-   **Breaking change**: The template parameters of
    [`CGAL::Tetrahedral_remeshing::Remeshing_vertex_base_3`](https://doc.cgal.org/5.6/Tetrahedral_remeshing/group__PkgTetrahedralRemeshingClasses.html#ga7ef4f8c0c1ed715c34389ea4ee851a92)
    and
    [`CGAL::Tetrahedral_remeshing::Remeshing_cell_base_3`](https://doc.cgal.org/5.6/Tetrahedral_remeshing/classCGAL_1_1Tetrahedral__remeshing_1_1Remeshing__cell__base__3.html)
    have been modified.

### [3D Mesh Generation](https://doc.cgal.org/5.6/Manual/packages.html#PkgMesh3)

-   Added two new named parameters to the named constructor
    [`CGAL::create_labeled_image_mesh_domain()`](https://doc.cgal.org/5.6/Mesh_3/classCGAL_1_1Labeled__mesh__domain__3.html#aec3f58e9883a8036a1b3e379df7d8fa9)
    for automatic detection and protection of 1D-curves that lie at the intersection of
    three or more subdomains extracted from labeled images.
-   Added
    [`CGAL::Sizing_field_with_aabb_tree`](https://doc.cgal.org/5.6/Mesh_3/structCGAL_1_1Sizing__field__with__aabb__tree.html),
    a geometry-aware sizing field for feature edges in polyhedral domains.
-   Added new meshing criterion
    [`edge_min_size`](https://doc.cgal.org/5.6/Mesh_3/classCGAL_1_1Mesh__criteria__3.html#a5f1c2649cb7ea346a3b6a2a8724b4df1)
    to avoid subdividing sharp edges that are shorter than a prescribed size bound.
-   Added new meshing criteria
    [`facet_min_size`](https://doc.cgal.org/5.6/Mesh_3/classCGAL_1_1Mesh__criteria__3.html#a5f1c2649cb7ea346a3b6a2a8724b4df1)
    and
    [`cell_min_size`](https://doc.cgal.org/5.6/Mesh_3/classCGAL_1_1Mesh__criteria__3.html#a5f1c2649cb7ea346a3b6a2a8724b4df1)
    to prevent Delaunay refinement from creating simplices smaller than a prescribed bound.
-   Deprecated usage of boost parameters in favor of function named parameters.

### [3D Periodic Mesh Generation](https://doc.cgal.org/5.6/Manual/packages.html#PkgPeriodic3Mesh3)

-   Periodic Mesh Generation now supports non-cubic domains.
-   Deprecated usage of boost parameters in favor of function named parameters.

### [Surface Mesh Simplification](https://doc.cgal.org/5.6/Manual/packages.html#PkgSurfaceMeshSimplification)
-   The stop predicates
    [`Count_stop_predicate`](https://doc.cgal.org/5.6/Surface_mesh_simplification/classCGAL_1_1Surface__mesh__simplification_1_1Count__stop__predicate.html)
    and
    [`Count_ratio_stop_predicate`](https://doc.cgal.org/5.6/Surface_mesh_simplification/classCGAL_1_1Surface__mesh__simplification_1_1Count__ratio__stop__predicate.html)
    are renamed to
    [`Edge_count_stop_predicate`](https://doc.cgal.org/5.6/Surface_mesh_simplification/classCGAL_1_1Surface__mesh__simplification_1_1Edge__count__stop__predicate.html)
    and
    [`Edge_count_ratio_stop_predicate`](https://doc.cgal.org/5.6/Surface_mesh_simplification/classCGAL_1_1Surface__mesh__simplification_1_1Edge__count__ratio__stop__predicate.html).
    Older versions have been deprecated.
-   Introduced
    [`Face_count_stop_predicate`](https://doc.cgal.org/5.6/Surface_mesh_simplification/classCGAL_1_1Surface__mesh__simplification_1_1Face__count__stop__predicate.html)
    and
    [`Face_count_ratio_stop_predicate`](https://doc.cgal.org/5.6/Surface_mesh_simplification/classCGAL_1_1Surface__mesh__simplification_1_1Face__count__ratio__stop__predicate.html),
    which can be used to stop the simplification algorithm based on a desired number of faces in the output, or a ratio between input and output face numbers.

### [2D Regularized Boolean Set Operations](https://doc.cgal.org/5.6/Manual/packages.html#PkgBooleanSetOperations2)
-   Exposed all required member functions of the
    [`GeneralPolygonWithHoles_2`](https://doc.cgal.org/5.6/Polygon/classGeneralPolygonWithHoles__2.html)
    concept (e.g.,
    [`clear_outer_boundary()`](https://doc.cgal.org/5.6/Polygon/classGeneralPolygonWithHoles__2.html#a9f5f035047505a2ccab3e68770f51bc6),
    [`clear_holes()`](https://doc.cgal.org/5.6/Polygon/classGeneralPolygonWithHoles__2.html#a2a507be648f127ac605da8c670ea2580),
    and
    [`clear()`](https://doc.cgal.org/5.6/Polygon/classGeneralPolygonWithHoles__2.html#a2ca4d9b43cc9216c1b2cdb080a915944)
    ).

[Release 5.5](https://github.com/CGAL/cgal/releases/tag/v5.5)
-----------

Release date: June 2022

### [3D Alpha Wrapping](https://doc.cgal.org/5.5/Manual/packages.html#PkgAlphaWrap3) (new package)

-   This component takes a 3D triangle mesh, soup, or point set as input, and generates a valid
    (watertight, intersection-free, and combinatorially 2-manifold) surface triangle mesh
    that contains the input.
    The algorithm proceeds by shrink-wrapping and refining a 3D Delaunay triangulation,
    starting from a loose bounding box of the input.
    Two user-defined parameters, alpha and offset, offer control over the maximum size of cavities
    where the shrink-wrapping process can enter, and the tightness of the final surface mesh
    to the input, respectively. Once combined, these parameters provide a means to trade fidelity
    to the input for complexity of the output.

    See also the [announcement page](https://www.cgal.org/2022/05/18/alpha_wrap/).

### [2D Straight Skeleton and Polygon Offsetting](https://doc.cgal.org/5.5/Manual/packages.html#PkgStraightSkeleton2) (breaking change)
-   Fix the output of the function [CGAL::create_exterior_skeleton_and_offset_polygons_with_holes_2()](https://doc.cgal.org/5.5/Straight_skeleton_2/group__PkgStraightSkeleton2OffsetFunctions.html#gaa159f093e5d6d7fdb62c1660a44f95fe)
    to not take into account the offset of the outer frame.
-   Fix the computation of the exterior offset of a polygon with holes that was not computing the offset of the holes

### [3D Convex Hulls](https://doc.cgal.org/5.5/Manual/packages.html#PkgConvexHull3)

-   Added an [overload of the function `CGAL::convex_hull_3()`](https://doc.cgal.org/5.5/Convex_hull_3/group__PkgConvexHull3Functions.html#ga52fca4745c2ef0351063fbe66b035fd1), which writes the result in an indexed triangle set.

### [2D Polygons](https://doc.cgal.org/5.5/Manual/packages.html#PkgPolygon2)

-   Add vertex, edge, and hole ranges.
-   The concept [`GeneralPolygonWithHoles_2`](https://doc.cgal.org/5.5/Polygon/classGeneralPolygonWithHoles__2.html) now requires the nested type `Polygon_2` instead of `General_polygon_2`.

### [2D Regularized Boolean Set-Operations](https://doc.cgal.org/5.5/Manual/packages.html#PkgBooleanSetOperations2)
-   The concept [`GeneralPolygonSetTraits_2`](https://doc.cgal.org/5.5/Boolean_set_operations_2/classGeneralPolygonSetTraits__2.html) now requires the nested type `Construct_polygon_with_holes_2` instead of `Construct_general_polygon_with_holes_2`.

### [Combinatorial Maps](https://doc.cgal.org/5.5/Manual/packages.html#PkgCombinatorialMaps)

-   Removed old code deprecated in CGAL 4.9 and 4.10 (global functions, and information associated with darts).

### [2D Arrangements](https://doc.cgal.org/5.5/Manual/packages.html#PkgArrangementOnSurface2)
-   Fixed the `intersect_2`, `compare_y_at_x_right`, and `compare_y_at_x_left` function objects of the traits class template [`Arr_geodesic_arc_on_sphere_traits_2`](https://doc.cgal.org/5.5/Arrangement_on_surface_2/classCGAL_1_1Arr__geodesic__arc__on__sphere__traits__2.html) that handles geodesic arcs on sphere and applied a small syntactical fix to the tracing traits.

### [Tetrahedral Mesh Generation](https://doc.cgal.org/5.5/Manual/packages.html#PkgMesh3)

-   Added the function
    [`remove_isolated_vertices()`](https://doc.cgal.org/5.5/Mesh_3/classCGAL_1_1Mesh__complex__3__in__triangulation__3.html#ace57c4e777da457c6e33b4f6e89949ce)
    as a post-processing step for the tetrahedral mesh generation.

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

[Release 5.4](https://github.com/CGAL/cgal/releases/tag/v5.4)
-----------

Release date: January 2022

### [General changes](https://doc.cgal.org/5.4/Manual/general_intro.html)

-   Added the cmake target `CGAL::CGAL_Basic_viewer` to ease the compilation of programs using
    the basic viewer-based function `CGAL::draw()`. This target will define the macro and link with
    `CGAL_Qt5` target when linked with it.

-   The kernel providing exact constructions and exact predicates
    ([`CGAL::Exact_predicates_exact_constructions_kernel`](https://doc.cgal.org/5.4/Kernel_23/classCGAL_1_1Exact__predicates__exact__constructions__kernel.html))
    is now thread-safe.
    See changes in `2D and 3D Linear Geometry Kernel` for more details.

-   The class `Geomview_stream` and all the dependent features have
    been removed from CGAL. Those features were actually no longer
    supported since CGAL-5.3 but it was not properly announced.

### [Shape Regularization](https://doc.cgal.org/5.4/Manual/packages.html#PkgShapeRegularization) (new package)

-   This package enables to regularize a set of segments and open or closed contours in 2D
    and a set of planes in 3D such that all input objects are rotated and aligned with respect to the
    user-specified conditions. In addition, it provides a global regularization framework that can be
    adjusted for the user needs and any type of geometric objects.

### [Weights](https://doc.cgal.org/5.4/Manual/packages.html#PkgWeights) (new package)

-   This package provides a simple and unified interface to different types of weights.
    In particular, it groups all weights into three category: analytic weights including
    all basic weights which can be computed analytically for a query point with respect to its
    local neighbors in 2D and 3D; barycentric weights, including all weights which can be computed
    for a query point with respect to the vertices of a planar polygon; and weighting regions,
    including all weights which are used to balance other weights.

### [2D Generalized Barycentric Coordinates](https://doc.cgal.org/5.4/Manual/packages.html#PkgBarycentricCoordinates2) (major changes)

-   **Breaking change**: The headers `Segment_coordinates_2.h` and `Triangle_coordinates_2.h` are
    renamed to `segment_coordinates_2.h` and `triangle_coordinates_2.h`.
-   The classes [`Segment_coordinates_2`](https://doc.cgal.org/5.4/Barycentric_coordinates_2/classCGAL_1_1Barycentric__coordinates_1_1Segment__coordinates__2.html)
    and [`Triangle_coordinates_2`](https://doc.cgal.org/5.4/Barycentric_coordinates_2/classCGAL_1_1Barycentric__coordinates_1_1Triangle__coordinates__2.html)
    are deprecated. The free functions [`compute_segment_coordinates_2()`](https://doc.cgal.org/5.4/Barycentric_coordinates_2/classCGAL_1_1Barycentric__coordinates_1_1Segment__coordinates__2.html#a134d363dccaeecb5621fa608fac76eaf)
    and [`compute_triangle_coordinates_2()`](https://doc.cgal.org/5.4/Barycentric_coordinates_2/classCGAL_1_1Barycentric__coordinates_1_1Triangle__coordinates__2.html#a958fee3ad9613d7bfa9d7a976aa3548f)
    are deprecated as well. Instead, the free functions [`segment_coordinates_2()`](https://doc.cgal.org/5.4/Barycentric_coordinates_2/group__PkgBarycentricCoordinates2RefFunctions.html#gab856ca68d37f58e6cdf74c8aac6f4245)
    and [`triangle_coordinates_2()`](https://doc.cgal.org/5.4/Barycentric_coordinates_2/group__PkgBarycentricCoordinates2RefFunctions.html#gaa378786f8996dbcefe7923ebb711e4dd)
    should be used.
-   The enums [`Query_point_location`](https://doc.cgal.org/5.4/Barycentric_coordinates_2/namespaceCGAL_1_1Barycentric__coordinates.html#aedeeb072a2024053a016afd15e591331)
    and [`Type_of_algorithm`](https://doc.cgal.org/5.4/Barycentric_coordinates_2/namespaceCGAL_1_1Barycentric__coordinates.html#a5e5682512438422f23d6080edc49c05b)
    are deprecated. Instead, the enum [`Computation_policy_2`](https://doc.cgal.org/5.4/Barycentric_coordinates_2/namespaceCGAL_1_1Barycentric__coordinates.html#a478bbcec416216b2274ee4b4e97b0e6c)
    should be used.
-   The classes [`Wachspress_2`](https://doc.cgal.org/5.4/Barycentric_coordinates_2/classCGAL_1_1Barycentric__coordinates_1_1Wachspress__2.html),
    [`Discrete_harmonic_2`](https://doc.cgal.org/5.4/Barycentric_coordinates_2/classCGAL_1_1Barycentric__coordinates_1_1Discrete__harmonic__2.html),
    [`Mean_value_2`](https://doc.cgal.org/5.4/Barycentric_coordinates_2/classCGAL_1_1Barycentric__coordinates_1_1Mean__value__2.html),
    and [`Generalized_barycentric_coordinates_2`](https://doc.cgal.org/5.4/Barycentric_coordinates_2/classCGAL_1_1Barycentric__coordinates_1_1Generalized__barycentric__coordinates__2.html)
    are deprecated. As consequence, the concept [`BarycentricCoordinates_2`](https://doc.cgal.org/5.4/Barycentric_coordinates_2/classCGAL_1_1Barycentric__coordinates_1_1BarycentricCoordinates__2.html)
    is deprecated as well. Instead, the classes [`Wachspress_coordinates_2`](https://doc.cgal.org/5.4/Barycentric_coordinates_2/classCGAL_1_1Barycentric__coordinates_1_1Wachspress__coordinates__2.html),
    [`Discrete_harmonic_coordinates_2`](https://doc.cgal.org/5.4/Barycentric_coordinates_2/classCGAL_1_1Barycentric__coordinates_1_1Discrete__harmonic__coordinates__2.html),
    and [`Mean_value_coordinates_2`](https://doc.cgal.org/5.4/Barycentric_coordinates_2/classCGAL_1_1Barycentric__coordinates_1_1Mean__value__coordinates__2.html)
    should be used.
-   Added the class [`Harmonic_coordinates_2`](https://doc.cgal.org/5.4/Barycentric_coordinates_2/classCGAL_1_1Barycentric__coordinates_1_1Harmonic__coordinates__2.html)
    to compute approximate harmonic coordinates in 2D.
    These coordinates satisfy all properties of barycentric coordinates inside any simple polygon.
-   Added a new concept [`DiscretizedDomain_2`](https://doc.cgal.org/5.4/Barycentric_coordinates_2/classCGAL_1_1Barycentric__coordinates_1_1DiscretizedDomain__2.html)
    and a model of this concept called [`Delaunay_domain_2`](https://doc.cgal.org/5.4/Barycentric_coordinates_2/classCGAL_1_1Barycentric__coordinates_1_1Delaunay__domain__2.html),
    which is based on the [Mesh 2](https://doc.cgal.org/5.4/Manual/packages.html#PkgMesh2) package.
    A model of this concept is required to use [`Harmonic_coordinates_2`](https://doc.cgal.org/5.4/Barycentric_coordinates_2/classCGAL_1_1Barycentric__coordinates_1_1Harmonic__coordinates__2.html).
-   Added free functions to compute Wachspress, discrete harmonic, and mean value coordinates.
-   All free functions and classes are now using ranges and property maps.

### [2D and 3D Linear Geometry Kernel](https://doc.cgal.org/5.4/Manual/packages.html#PkgKernel23)

-   Most operations on [`CGAL::Exact_predicates_exact_constructions_kernel`](https://doc.cgal.org/5.4/Kernel_23/classCGAL_1_1Exact__predicates__exact__constructions__kernel.html)
    objects are now thread-safe if [`CGAL::Exact_rational`](https://doc.cgal.org/5.4/Number_types/group__nt__cgal.html#ga0849ff44771b19582218ebdfa5614f64)
    is [`mpq_class`](https://doc.cgal.org/5.3/Number_types/classmpq__class.html) (from `GMPXX`),
    `boost::multiprecision::mpq_rational`
    or [`CGAL::Quotient<CGAL::MP_Float>`](https://doc.cgal.org/5.3/Number_types/classCGAL_1_1MP__Float.html).
    The objects are not atomic though, so the usual restrictions on avoiding race conditions apply.
    For users who do not use threads, this can be disabled with `CGAL_HAS_NO_THREADS`.

-   Added documentation for the class [`Projection_traits_3`](https://doc.cgal.org/5.4/Kernel_23/classCGAL_1_1Projection__traits__3.html),
    which enables the use of 2D algorithms on the projections of 3D data onto an arbitrary plane.

-   Added `construct_centroid_2_object()` and `compute_determinant_2_object()`
    in [`Projection_traits_xy_3`](https://doc.cgal.org/5.4/Kernel_23/classCGAL_1_1Projection__traits__xy__3.html),
    [`Projection_traits_xz_3`](https://doc.cgal.org/5.4/Kernel_23/classCGAL_1_1Projection__traits__xz__3.html),
    and [`Projection_traits_yz_3`](https://doc.cgal.org/5.4/Kernel_23/classCGAL_1_1Projection__traits__yz__3.html)
    classes.

-   Added the functor
    [`NonZeroCoordinateIndex_3`](https://doc.cgal.org/5.4/Kernel_23/classKernel_1_1NonZeroCoordinateIndex__3.html)
    to the concept [`Kernel`](https://doc.cgal.org/5.4/Kernel_23/classKernel.html) with `int operator()(Vector_3)`
    which returns the index of any coordinate of the vector different from zero, or `-1`.

### [dD Kernel](https://doc.cgal.org/5.4/Manual/packages.html#PkgKernelD)

-   Most operations on [`Epeck_d`](https://doc.cgal.org/5.4/Kernel_d/structCGAL_1_1Epeck__d.html)
    objects are now thread-safe, see 2D and 3D Linear Geometry Kernel for details.

### [2D Arrangements](https://doc.cgal.org/5.4/Manual/packages.html#PkgArrangementOnSurface2)

-   **Breaking Change:** The traits function objects `Compare_x_at_limit_2` and `Compare_x_near_limit_2`
    are renamed to `Compare_x_on_boundary_2` and `Compare_x_near_boundary_2`, respectively.

-   A [new hierarchy of traits concepts](https://doc.cgal.org/5.4/Arrangement_on_surface_2/group__PkgArrangementOnSurface2Concepts.html)
    has been introduced.
    It captures all the valid combinations of boundary conditions for the 4 boundary sides of the parameter space.
    The 4 boundaries are Bottom, Top, Left, and Right. Each boundary side can be either contracted, identified, close, open, or oblivious.
    Not all possible combinations are valid. If one side is identified then the other must be as well. Two adjacent sides cannot be contracted.

-   A new geometric traits, [`Arr_geodesic_arc_on_sphere_traits_2`](https://doc.cgal.org/5.4/Arrangement_on_surface_2/classCGAL_1_1Arr__geodesic__arc__on__sphere__traits__2.html)
    has been introduced. It handles arcs of great circles embedded on the unit sphere.

### [2D Regularized Boolean Set-Operations](https://doc.cgal.org/5.4/Manual/packages.html#PkgBooleanSetOperations2)

-   Added an extra parameter (`UsePolylines`) to all free functions (
    [`complement()`](https://doc.cgal.org/5.4/Boolean_set_operations_2/group__boolean__complement.html),
    [`do_intersect()`](https://doc.cgal.org/5.4/Boolean_set_operations_2/group__boolean__do__intersect.html),
    [`intersection()`](https://doc.cgal.org/5.4/Boolean_set_operations_2/group__boolean__intersection.html),
    [`join()`](https://doc.cgal.org/5.4/Boolean_set_operations_2/group__boolean__join.html),
    [`difference()`](https://doc.cgal.org/5.4/Boolean_set_operations_2/group__boolean__difference.html),
    [`symmetric_difference()`](https://doc.cgal.org/5.4/Boolean_set_operations_2/group__boolean__symmetric__difference.html),
    and [`oriented_side`](https://doc.cgal.org/5.4/Boolean_set_operations_2/group__boolean__oriented__side.html))
    to control whether to use `Arr_polyline_traits_2` as default traits. It is the new default as it provides better performances in general.

### [3D Mesh Generation](https://doc.cgal.org/5.4/Manual/packages.html#PkgMesh3)

-   Added support of weighted images for an improved quality of meshes generated from labeled images,
    along with a function [`CGAL::Mesh_3::generate_label_weights()`](https://doc.cgal.org/5.4/Mesh_3/namespaceCGAL_1_1Mesh__3.html#ae5914bf77180ff8948c08046154ee727)
    to generate the weights.

### [Polygon Mesh Processing](https://doc.cgal.org/5.4/Manual/packages.html#PkgPolygonMeshProcessing)

-   Added the function [`CGAL::Polygon_mesh_processing::match_faces()`](https://doc.cgal.org/5.4/Polygon_mesh_processing/group__measure__grp.html#ga10f7cd81645bafe936ac5eb4e58e67ef),
    which, given two polygon meshes, identifies their common faces as well as faces present in only either of them.

-   Added the functions: [`CGAL::Polygon_mesh_processing::bounded_error_Hausdorff_distance()`](https://doc.cgal.org/5.4/Polygon_mesh_processing/group__PMP__distance__grp.html#ga6d4ecea831c33ac10eec42b5021fc183)
    that computes an estimate of the one-sided Hausdorff distance between two triangle meshes which
    is bounded by a user-specified error bound; [`CGAL::Polygon_mesh_processing::bounded_error_symmetric_Hausdorff_distance()`](https://doc.cgal.org/5.4/Polygon_mesh_processing/group__PMP__distance__grp.html#ga9a7a682b5d9523135c8502e72117dffd)
    that computes an estimate of the symmetric Hausdorff distance bounded by a user-specified error bound;
    and [`CGAL::Polygon_mesh_processing::is_Hausdorff_distance_larger()`](https://doc.cgal.org/5.4/Polygon_mesh_processing/group__PMP__distance__grp.html#gab19e751107025a443e86baa9763aebf3)
    that returns `true` if the bounded-error Hausdorff distance between two meshes is larger than the user-specified
    max distance.

-   Added the functions [`CGAL::Polygon_mesh_processing::squared_edge_length()`](https://doc.cgal.org/5.4/Polygon_mesh_processing/group__measure__grp.html#ga30fa03722cd7aa599f6dcb115f54fec5)
    and [`CGAL::Polygon_mesh_processing::squared_face_area()`](https://doc.cgal.org/5.4/Polygon_mesh_processing/group__measure__grp.html#ga6eda3738815fd678df225f79ccfc3e03),
    which, compared to [`CGAL::Polygon_mesh_processing::edge_length()`](https://doc.cgal.org/5.4/Polygon_mesh_processing/group__measure__grp.html#gae1674775d9fecada7f25710f425cff3a)
    and [`CGAL::Polygon_mesh_processing::face_area()`](https://doc.cgal.org/5.4/Polygon_mesh_processing/group__measure__grp.html#ga6a1d7a825c09490b1e6613295343482b),
    enable avoiding square-root operations.

-   Added more functions in the [visitor of the corefinement based methods](https://doc.cgal.org/5.4/Polygon_mesh_processing/classPMPCorefinementVisitor.html)
    to track all vertex creations.

-   Added an option to [`CGAL::Polygon_mesh_processing::self_intersections()`](https://doc.cgal.org/5.4/Polygon_mesh_processing/group__PMP__intersection__grp.html#gaf19c80ec12cbff7ebe9e69453f1d40b8)
    to report only a limited number of intersections (`maximum_number()`).

### [The Heat Method](https://doc.cgal.org/5.4/Manual/packages.html#PkgHeatMethod)

-   **Breaking change**: Added the functor `Compute_squared_length_3` providing `operator(const Vector_3& v)`,
    which computes the squared length of `v`, to the [`HeatMethodTraits_3`](https://doc.cgal.org/5.4/Heat_method_3/classHeatMethodTraits__3.html)
    concept.

### [Point Set Processing](https://doc.cgal.org/5.4/Manual/packages.html#PkgPointSetProcessing3)

-   Added support for [`libpointmatcher::GenericDescriptorOutlierFilter`](https://github.com/ethz-asl/libpointmatcher)
    that enables providing a map from a point to a weight associated with this point.

### [Shape Detection](https://doc.cgal.org/5.4/Manual/packages.html#PkgShapeDetection)

-   Added new shapes to the Region Growing algorithm on a point set: circles in 2D, spheres in 3D,
    and cylinders in 3D.

###  [CGAL and Solvers](https://doc.cgal.org/5.4/Manual/packages.html#PkgSolverInterface)

-   Added support for the [OSQP solver](https://osqp.org/). This solver enables to efficiently compute
    the convex Quadratic Programming (QP) problems arising in the context of several packages.


[Release 5.3](https://github.com/CGAL/cgal/releases/tag/v5.3)
-----------

Release date: July 2021

### [General changes](https://doc.cgal.org/5.3/Manual/general_intro.html)

-   The support for the compiled version of CGAL is dropped. Only the header-only version is supported.

-   On Windows, the type used for [`CGAL::Exact_rational`](https://doc.cgal.org/5.3/Number_types/group__nt__cgal.html#ga0849ff44771b19582218ebdfa5614f64),
    in `Epick` and indirectly (through [`Lazy_exact_nt`](https://doc.cgal.org/5.3/Number_types/classCGAL_1_1Lazy__exact__nt.html))
   `Epeck` may now be `boost::multiprecision::mpq_rational`, as has been the case on other platforms
   for several releases. This depends on various options and is added to a list that includes
   [`mpq_class`](https://doc.cgal.org/5.3/Number_types/classmpq__class.html),
   [`CGAL::Gmpq`](https://doc.cgal.org/5.3/Number_types/classCGAL_1_1Gmpq.html),
   [`leda_rational`](https://doc.cgal.org/5.3/Number_types/classleda__rational.html)
   and [`CGAL::Quotient<CGAL::MP_Float>`](https://doc.cgal.org/5.3/Number_types/classCGAL_1_1MP__Float.html).

### [Quadtrees, Octrees, and Orthtrees](https://doc.cgal.org/5.3/Manual/packages.html#PkgOrthtree) (new package)

-   This package implements a tree data structure in which each node encloses a hypercubic section
    of space and each non-leave node has hypercubic children whose edge lengths are half its edge length.
    Such a data structure is known as a quadtree in 2D, an octree in 3D, and is generalized
    as an "orthtree" in higher dimensions.

### [Triangulations on the Sphere](https://doc.cgal.org/5.3/Manual/packages.html#PkgTriangulationOnSphere2) (new package)

-   This package enables the construction and manipulation of Delaunay triangulations on the 2-sphere.
    Triangulations are built incrementally and can be modified by insertion or removal of vertices.
    Point location querying and primitives to build the dual Voronoi diagram are provided.

### File Input / Output

-   Point set, polygon soup, and polygon mesh file I/O functions have been harmonized and documented:
    -   Point set I/O functions can be found in the packages [Point_set_processing_3](https://doc.cgal.org/5.3/Manual/packages.html#PkgPolygonMeshProcessing), and [Point_set_3](https://doc.cgal.org/5.3/Manual/packages.html#PkgPointSet3).
    -   Polygon mesh I/O functions can be found in the package [BGL](https://doc.cgal.org/5.3/Manual/packages.html#PkgBGL).
    -   Polygon soup I/O can be found in the package [Stream_support](https://doc.cgal.org/5.3/Manual/packages.html#PkgStreamSupport).

A comprehensive list of the supported file formats is available in the Stream_support package
[here](https://doc.cgal.org/5.3/Stream_support/index.html#IOstreamSupportedFormats);
inversely, the following [page](https://doc.cgal.org/5.3/Stream_support/IOStreamSupportedFileFormats.html)
can be used to find out which CGAL data structures can be used given a specific file format.

### [Requirements](https://doc.cgal.org/5.3/Manual/thirdparty.html)

-   The CMake minimal version is now `3.14`.
-   The GNU compiler g++ versions 6 and 7 are no longer tested. Only version 8.3 or later are supported

### [2D and 3D Linear Geometry Kernel](https://doc.cgal.org/5.3/Manual/packages.html#PkgKernel23)

-   Added `is_translation()`, `is_scaling()`, `is_reflection()`, and `is_rotation()` to the classes
    [`Aff_transformation_2`](https://doc.cgal.org/5.3/Kernel_23/classCGAL_1_1Aff__transformation__2.html)
    and [`Aff_transformation_3`](https://doc.cgal.org/5.3/Kernel_23/classCGAL_1_1Aff__transformation__3.html),
    which enable determining if the transformations use a specialized representation internally.

### [2D Regularized Boolean Set-Operations](https://doc.cgal.org/5.3/Manual/packages.html#PkgBooleanSetOperations2)
-   Added documentation for the free functions [`oriented_side(const Point_2& p, ....)`](https://doc.cgal.org/5.3/Boolean_set_operations_2/group__boolean__oriented__side.html)
    that accept a point and a polygon.
-   Documentation has been improved across the whole package.

### [Polygon Mesh Processing](https://doc.cgal.org/5.3/Manual/packages.html#PkgPolygonMeshProcessing)

-   Added the class [`CGAL::Polyhedral_envelope`](https://doc.cgal.org/5.3/Polygon_mesh_processing/structCGAL_1_1Polyhedral__envelope.html),
    providing a way to quickly check if a primitive (point, segment, or triangle)
    is within a polyhedral envelope around a set of triangles. It is based on the work of
    Bolun Wang, Teseo Schneider, Yixin Hu, Marco Attene, and Daniele Panozzo.
    "Exact and efficient polyhedral envelope containment check." (ACM Trans. Graph., 39-4, July 2020).
-   Added more functions in the [visitor of the corefinement based methods](https://doc.cgal.org/5.3/Polygon_mesh_processing/classPMPCorefinementVisitor.html)
    to track all edge creations.

### [Surface Mesh Topology](https://doc.cgal.org/5.3/Manual/packages.html#PkgSurfaceMeshTopologySummary)
-   Added the function [`CGAL::Surface_mesh_topology::Curves_on_surface_topology::is_homotopic_to_simple_cycle()`](https://doc.cgal.org/5.3/Surface_mesh_topology/classCGAL_1_1Surface__mesh__topology_1_1Curves__on__surface__topology.html#a8d7c4cba2cf2cff542f5cd93117233db),
    which can be used to determine whether a closed path on a surface mesh can be continuously
    transformed to a cycle without self intersection.

### [Surface Mesh Simplification](https://doc.cgal.org/5.3/Manual/packages.html#PkgSurfaceMeshSimplification)
-   Added a filtering mechanism so that costly tests get only applied to the next candidate for the edge collapse.
-   Added the class [`Polyhedral_envelope_filter`](https://doc.cgal.org/5.3/Surface_mesh_simplification/classCGAL_1_1Surface__mesh__simplification_1_1Polyhedral__envelope__filter.html),
    which enables to perform mesh simplification  inside a polyhedral envelope of the input mesh.

### [2D Polyline Simplification](https://doc.cgal.org/5.3/Manual/packages.html#PkgPolylineSimplification2)
-   When polylines have common subsequences of vertices, these subsequences may now be simplifified simultaneously.

### [dD Triangulations](https://doc.cgal.org/5.3/Manual/packages.html#PkgTriangulations)
-   Added the function [`insert_if_in_star()`](https://doc.cgal.org/5.3/Triangulation/classCGAL_1_1Regular__triangulation.html#aa8df2d138f341939e834bcdd7cb6c71a)
    to the class [`CGAL::Regular_triangulation`](https://doc.cgal.org/5.3/Triangulation/classCGAL_1_1Regular__triangulation.html),
    which enables users to insert a point `p` in a regular triangulation on the condition that `p`
    appears post-insertion in the star of a user-specified, existing vertex.

### [2D and 3D Alpha Shapes](https://doc.cgal.org/5.3/Manual/packages.html#PkgAlphaShapes2)
-   **Breaking change**: The following deprecated classes have been removed: `Alpha_shape_euclidean_traits_2`,
    `Weighted_alpha_shape_euclidean_traits_2`, `Alpha_shape_euclidean_traits_3`, and
    `Weighted_alpha_shape_euclidean_traits_3`. All CGAL kernel can be used directly as models
    of the concepts of the 2D and 3D Alpha Shape packages.

### [Classification](https://doc.cgal.org/5.3/Manual/packages.html#PkgClassification)
-   **Breaking change**: the support for TensorFlow has been dropped; the
    classifier `CGAL::TensorFlow::Neural_network_classifier` has been removed.


[Release 5.2](https://github.com/CGAL/cgal/releases/tag/v5.2)
-----------

Release date: December 2020

### [dD Geometry Kernel](https://doc.cgal.org/5.2/Manual/packages.html#PkgKernelD)

-   The kernels [`Epick_d`](https://doc.cgal.org/5.2/Kernel_d/structCGAL_1_1Epick__d.html)
    and [`Epeck_d`](https://doc.cgal.org/5.2/Kernel_d/structCGAL_1_1Epeck__d.html) gain two new functors:
    [`Compute_power_product_d`](https://doc.cgal.org/5.2/Kernel_d/classCGAL_1_1Epeck__d_1_1Compute__power__product__d.html)
    and [`Construct_power_sphere_d`](https://doc.cgal.org/5.2/Kernel_d/classCGAL_1_1Epeck__d_1_1Construct__power__sphere__d.html),
    to deal with weighted points.

### [CGAL and the Boost Graph Library (BGL)](https://doc.cgal.org/5.2/Manual/packages.html#PkgBGL)

-   Added a convenience header, [`CGAL/boost/graph/graph_traits_inheritance_macros.h`](https://doc.cgal.org/5.2/BGL/graph__traits__inheritance__macros_8h.html),
    which enables easily making any class inheriting from a model of a face graph concept, a model of the same concept.
-   Added the function [`can_add_face()`](https://doc.cgal.org/5.2/BGL/group__PkgBGLEulerOperations.html#ga7dc63595108097b6e28b04fe962135f0),
    which tests whether a new face defined by a range of vertices can be added.

### [3D Fast Intersection and Distance Computation (AABB Tree)](https://doc.cgal.org/5.2/Manual/packages.html#PkgAABBTree)

-   Added the move constructor and the assignment operator to the
    [AABB Tree](https://doc.cgal.org/5.2/AABB_tree/classCGAL_1_1AABB__tree.html) class.

### [2D Arrangements](https://doc.cgal.org/5.2/Manual/packages.html#PkgArrangementOnSurface2)

-   Replaced the use of legacy
    [`CGAL::Object`](https://doc.cgal.org/5.2/STL_Extension/classCGAL_1_1Object.html)
    to modern `boost::variant`.
-   Changed make-x-monotone return type from legacy
    [`CGAL::Object`](https://doc.cgal.org/5.2/STL_Extension/classCGAL_1_1Object.html)
    to `boost::variant` in all traits concepts and models.
    As there exists an implicit conversion from `boost::variant` to `CGAL::Object`,
    the new code is backward compatible. However, it is recommended that all calls
    to the make-x-monotone functions are fixed to use the new return type.
-   Changed [`decompose()`](https://doc.cgal.org/5.2/Arrangement_on_surface_2/group__PkgArrangementOnSurface2Funcs.html#gae20b2917f6de15db9bf025f83abf8e89)
    interface to use `boost::variant` instead of legacy
    [`CGAL::Object`](https://doc.cgal.org/5.2/STL_Extension/classCGAL_1_1Object.html)
    As explained above, the code is backward compatible. However, it is recommended
    that all calls to `decompose()` are fixed to use the new interface.

### [Surface Mesh](https://doc.cgal.org/5.2/Manual/packages.html#PkgSurfaceMesh)

-   Added the function [`clear_without_removing_property_maps()`](https://doc.cgal.org/5.2/Surface_mesh/classCGAL_1_1Surface__mesh.html#aad000a07a5ada30536f194b28b59d111)
    to clear a mesh but keep all the created property maps added.
-   Added the functions [`remove_property_maps<Index_type>()`](https://doc.cgal.org/5.2/Surface_mesh/classCGAL_1_1Surface__mesh.html#a2a3dd8c01f7fba7b640d85bfd1c41d90)
    and [`remove_all_property_maps()`](https://doc.cgal.org/5.2/Surface_mesh/classCGAL_1_1Surface__mesh.html#a5696da09300f3d0eafed117668bb3bec)
    to remove all added property maps by index type or all of them respectively.
-   Added the functions [`set_recycle_garbage()`](https://doc.cgal.org/5.2/Surface_mesh/classCGAL_1_1Surface__mesh.html#a40ada5068bf6d529a511c46767dfd21d)
    and [`does_recycle_garbage()`](https://doc.cgal.org/5.2/Surface_mesh/classCGAL_1_1Surface__mesh.html#a081a87aaf7e56e6b4f9afba99967f8f4)
    to the class `Surface_mesh`.

### [Polygon Mesh Processing](https://doc.cgal.org/5.2/Manual/packages.html#PkgPolygonMeshProcessing)

-   Added a visitor to the functions
    [`CGAL::Polygon_mesh_processing::triangulate_face()`](https://doc.cgal.org/5.2/Polygon_mesh_processing/group__PMP__meshing__grp.html#ga70d65044f8c7309c24ade88fa280124a)
    and [`CGAL::Polygon_mesh_processing::triangulate_faces()`](https://doc.cgal.org/5.2/Polygon_mesh_processing/group__PMP__meshing__grp.html#gacaaff4d520500c530d9c3d5ebe2a0760),
    that enables the user to keep track of the newly created faces through the triangulation process.
-   Added an option in [`CGAL::Polygon_mesh_processing::corefine()`](https://doc.cgal.org/5.2/Polygon_mesh_processing/group__PMP__corefinement__grp.html#ga6447dee822aaf92016f34512ce0b3456),
    [`CGAL::Polygon_mesh_processing::split()`](https://doc.cgal.org/5.2/Polygon_mesh_processing/group__PMP__corefinement__grp.html#gaa491feee9e41f725332bea0ea1215578)
    and [`CGAL::Polygon_mesh_processing::clip()`](https://doc.cgal.org/5.2/Polygon_mesh_processing/group__PMP__corefinement__grp.html#ga30082762ba2d947cba304e2884d96a99)
    functions, which enable the operations to be performed on a mesh with
    self-intersections present in the intersection area.
-   Added an optional range parameter to [`CGAL::Polygon_mesh_processing::stitch_borders()`](https://doc.cgal.org/5.2/Polygon_mesh_processing/group__PMP__repairing__grp.html#ga8ae4352e67d2b099994ac8990c13bd41),
    which can be used to specify which boundary cycles are eligible for stitching.

### [Surface Mesh Parameterization](https://doc.cgal.org/5.2/Manual/packages.html#PkgSurfaceMeshParameterization)

-   Added a new parameterization method, [Iterative Authalic Parameterization](https://doc.cgal.org/5.2/Surface_mesh_parameterization/index.html#title11).
    It is based on the work of Jain, Hardik, Manuel Wollhaf, and Olaf Hellwich,
    "Learning to Reconstruct Symmetric Shapes using Planar Parameterization of 3D Surface."
    (IEEE International Conference on Computer Vision Workshops, 2019).

### [Classification](https://doc.cgal.org/5.2/Manual/packages.html#PkgClassification)

-   **Breaking change**: new IO format for the [`ETHZ::Random_Forest`](https://doc.cgal.org/5.2/Classification/classCGAL_1_1Classification_1_1ETHZ_1_1Random__forest__classifier.html) classifier:
    a conversion function from the outdated format to the new one is provided.

-   Added new functions to the class [`CGAL::Classification::Evaluation`](https://doc.cgal.org/5.2/Classification/classCGAL_1_1Classification_1_1Evaluation.html):
    [`append()`](https://doc.cgal.org/5.2/Classification/classCGAL_1_1Classification_1_1Evaluation.html#a20c5fc43af44c96ce0cae40375be934f)
    to enrich the evaluation with additional results;
    [`confusion()`](https://doc.cgal.org/5.2/Classification/classCGAL_1_1Classification_1_1Evaluation.html#a706a85bb1deefee350ce71855bc023e9)
    to access the confusion matrix;
    output functions to save the evaluation to and `ASCII` or `HTML` stream.
-   Added a new operator, [`CGAL::Classification::feature_cast<>`](https://doc.cgal.org/5.2/Classification/group__PkgClassificationFeature.html#gaf4b1504270f25061f63f05743a17e5d1),
    for easy conversions.
-   The classes [`CGAL::Classification::Feature_set`](https://doc.cgal.org/5.2/Classification/classCGAL_1_1Classification_1_1Feature__set.html)
    and [`CGAL::Classification::Label_set`](https://doc.cgal.org/5.2/Classification/classCGAL_1_1Classification_1_1Label__set.html)
    are now models of the concept [`Range`](https://doc.cgal.org/5.2/Circulator/classRange.html).
-   The class [`CGAL::Classification::Label`](https://doc.cgal.org/5.2/Classification/classCGAL_1_1Classification_1_1Label.html)
    now has attributes `index`, `standard_index` and `color`,
    with automatic selection if the ASPRS standard names are used.
-   Added new functions in [`CGAL::Classification::Point_set_feature_iterator`](https://doc.cgal.org/5.2/Classification/classCGAL_1_1Classification_1_1Point__set__feature__generator.html),
    to enable users to select which features should be generated.
-   Added a new function, [`CGAL::Classification::Label_set::is_valid_ground_truth()`](https://doc.cgal.org/5.2/Classification/classCGAL_1_1Classification_1_1Label__set.html#adeb3b046f640c091b1f123e982386e43),
    to help users check if a ground truth matches a given label set.

### [Point Set Processing](https://doc.cgal.org/5.2/Manual/packages.html#PkgPointSetProcessing3)

-   Added a function [`CGAL::scanline_orient_normals()`](https://doc.cgal.org/5.2/Point_set_processing_3/group__PkgPointSetProcessing3Algorithms.html#ga221d4efde44f42aefe153cb927138efe),
    which orients a point cloud by estimating a line of sight.

### [3D Convex Hulls](https://doc.cgal.org/5.2/Manual/packages.html#PkgConvexHull3)

-   Added the function [`CGAL::halfspace_intersection_interior_point_3()`](https://doc.cgal.org/5.2/Convex_hull_3/group__PkgConvexHull3Functions.html#ga9a1ead3126e42fbf46ef269466cddc8f),
    which can be used to retrieve the point that is the most interior a convex closed volume
    defined by the intersection of a set of halfspaces.

### [3D Triangulations](https://doc.cgal.org/5.2/Manual/packages.html#PkgTriangulation3)

-   Added new classes and functions to visit the cells and simplices intersected by a line segment,
    see Sections [Segment Cell Iterator](https://doc.cgal.org/5.2/Triangulation_3/classCGAL_1_1Triangulation__3.html#amgrp0d087ed77bb99ca595c92d2cd2ab59b9)
    and [Segment Simplex Iterator](https://doc.cgal.org/5.2/Triangulation_3/classCGAL_1_1Triangulation__3.html#amgrp2447c1d2dce281951a0a4d8aecd3f35d), respectively.


[Release 5.1](https://github.com/CGAL/cgal/releases/tag/v5.1)
-----------

Release date: September 2020

### [Tetrahedral Remeshing](https://doc.cgal.org/5.1/Manual/packages.html#PkgTetrahedralRemeshing) (new package)

-   This package implements a tetrahedral isotropic remeshing algorithm,
    that improves the quality of tetrahedra in terms of dihedral angles,
    while targeting a given edge length.

    See also the associated [blog entry](https://www.cgal.org/2020/08/07/Tetrahedral-remeshing/).

### [Surface Mesh Topology](https://doc.cgal.org/5.1/Manual/packages.html#PkgSurfaceMeshTopologySummary) (new package)

-   This package enables the computation of some topological invariants of surfaces, such as:
    - test if two (closed) curves on a combinatorial surface are homotopic. Users can choose
      between free homotopy and homotopy with fixed endpoints;
    - test is a curve is contractible;
    - compute shortest non-contractible cycles on a surface, with or without weights on edges.

    See also the associated [blog entry](https://www.cgal.org/2020/05/08/Surface_mesh_topology/).

### [Optimal Bounding Box](https://doc.cgal.org/5.1/Manual/packages.html#PkgOptimalBoundingBox) (new package)

-   This package implements an optimization algorithm that aims to construct a close approximation
    of the *optimal bounding box* of a mesh or a point set, which is defined as the smallest
    (in terms of volume) bounding box that contains a given mesh or point set.

    See also the associated [blog entry](https://www.cgal.org/2020/04/20/Optimal_bounding_box/).

### Installation

-   The CGAL\_Core library no longer requires `Boost.Thread`, even if the g++ compiler is used.
-   The minimal supported version of Boost is now 1.66.0.

### [Tutorials](https://doc.cgal.org/5.1/Manual/tutorials.html)

-   Two new, detailed tutorials have been added:
    - [Surface Reconstruction from Point Clouds](https://doc.cgal.org/5.1/Manual/tuto_reconstruction.html),
      which goes over a typical full processing pipeline in a CGAL environment.
    - [Geographic Information Systems (GIS)](https://doc.cgal.org/5.1/Manual/tuto_gis.html),
      which demonstrates usage of CGAL data structures and algorithms in the context of a typical GIS application.

    Both tutorials provide complete code.

### [2D and 3D Linear Geometry Kernel](https://doc.cgal.org/5.1/Manual/packages.html#PkgKernel23)

-   Added the functor [`CompareSignedDistanceToLine_2`](https://doc.cgal.org/5.1/Kernel_23/classKernel_1_1CompareSignedDistanceToLine__2.html)
    to the 2D/3D [`Kernel`](https://doc.cgal.org/5.1/Kernel_23/classKernel.html) concept to compare
    the signed distance of two points to a line, or the line passing through two given points.
    Corresponding functors in the model ([`Compare_signed_distance_to_line_2`](https://doc.cgal.org/5.1/Kernel_23/classKernel.html#a066d07dd592ac36ba7ee90988abd349f)) are also added.

### [dD Geometry Kernel](https://doc.cgal.org/5.1/Manual/packages.html#PkgKernelD)

-   The kernels [`Epick_d`](https://doc.cgal.org/5.1/Kernel_d/structCGAL_1_1Epick__d.html)
    and [`Epeck_d`](https://doc.cgal.org/5.1/Kernel_d/structCGAL_1_1Epeck__d.html) gain two new functors:
    [`Power_side_of_bounded_power_sphere_d`](https://doc.cgal.org/5.1/Kernel_d/classCGAL_1_1Epeck__d_1_1Power__side__of__bounded__power__sphere__d.html)
    and [`Compute_squared_radius_smallest_orthogonal_sphere_d`](https://doc.cgal.org/5.1/Kernel_d/classCGAL_1_1Epeck__d_1_1Compute__squared__radius__smallest__orthogonal__sphere__d.html).
    Those are essential for the computation of weighted alpha-complexes.

### [Surface Mesh](https://doc.cgal.org/5.1/Manual/packages.html#PkgSurfaceMesh)

-   **Breaking change**: The function [`CGAL::Surface_mesh::clear()`](https://doc.cgal.org/5.1/Surface_mesh/classCGAL_1_1Surface__mesh.html#a247d4ad3e6b106ae22e5306203812642)
    now removes all non-default properties instead of just emptying them.

### [CGAL and the Boost Graph Library (BGL)](https://doc.cgal.org/5.1/Manual/packages.html#PkgBGL)

-   Added the function [`CGAL::alpha_expansion_graphcut()`](https://doc.cgal.org/5.1/BGL/group__PkgBGLPartition.html#ga79c3f58b577af51d1140450729d38f22),
    which regularizes a multi-label partition over a user-defined graph.
-   Added the function [`CGAL::regularize_face_selection_borders()`](https://doc.cgal.org/5.1/BGL/group__PkgBGLSelectionFct.html#gac71322b0cc7d7d59447531d5e5e345b6),
    which uses this alpha expansion graphcut to regularize the borders of a selected faces on a triangle mesh.
-   Added the function [`CGAL::set_triangulation_ids()`](https://doc.cgal.org/5.1/BGL/group__BGLGraphExternalIndices.html#ga1a22cf8bdde32fcdf1a4a78966eed630),
    which must be used to initialize vertex, edge, and face indices of a triangulation meant to be used with BGL algorithms.

### [3D Fast Intersection and Distance Computation (AABB Tree)](https://doc.cgal.org/5.1/Manual/packages.html#PkgAABBTree)

-   The behavior of the internal search tree used to accelerate distance queries has changed:
    usage of the internal search tree will now be enabled by default, and its construction
    will be triggered by the first distance query. Automatic construction and usage can be disabled
    by calling [`CGAL::AABB_tree::do_not_accelerate_distance_queries()`](https://doc.cgal.org/5.1/AABB_tree/classCGAL_1_1AABB__tree.html#abde62f52ccdf411847151aa5000ba4a4)
    before the first distance query, and the tree can be built at any moment by calling
    [`CGAL::AABB_tree::accelerate_distance_queries()`](https://doc.cgal.org/5.1/AABB_tree/classCGAL_1_1AABB__tree.html#a5d3877d3f2afbd09341eb4b8c230080b).
-   **Breaking change**: [`CGAL::AABB_tree::accelerate_distance_queries()`](https://doc.cgal.org/5.1/AABB_tree/classCGAL_1_1AABB__tree.html#a5d3877d3f2afbd09341eb4b8c230080b)
    and [`CGAL::AABB_tree::do_not_accelerate_distance_queries()`](https://doc.cgal.org/5.1/AABB_tree/classCGAL_1_1AABB__tree.html#abde62f52ccdf411847151aa5000ba4a4)
    are no longer `const` functions.

### [2D Arrangements](https://doc.cgal.org/5.1/Manual/packages.html#PkgArrangementOnSurface2)

 -   Changed intersection return type from legacy [`CGAL::Object`](https://doc.cgal.org/5.1/STL_Extension/classCGAL_1_1Object.html)
     to modern `boost::variant` in all traits concepts and models.
     As there exists an implicit conversion from `boost::variant` to `CGAL::Object`, the
     new code is backward compatible. However, it is recommended that all calls
     to the intersection functions are fixed to use the new return type.

### [2D Regularized Boolean Set-Operations](https://doc.cgal.org/5.1/Manual/packages.html#PkgBooleanSetOperations2)

 -   Changed intersection return type from legacy [`CGAL::Object`](https://doc.cgal.org/5.1/STL_Extension/classCGAL_1_1Object.html)
     to modern `boost::variant` in the concept [`ArrDirectionalTraits::Intersect_2`](https://doc.cgal.org/5.1/Boolean_set_operations_2/namespaceArrDirectionalTraits.html)
     and its models.

### [2D Minkowski Sums](https://doc.cgal.org/5.1/Manual/packages.html#PkgMinkowskiSum2)

 -   Changed intersection return type from legacy [`CGAL::Object`](https://doc.cgal.org/5.1/STL_Extension/classCGAL_1_1Object.html)
     to modern `boost::variant` in the (internally used) model `Arr_labeled_traits_2`.

### [dD Spatial Searching](https://doc.cgal.org/5.1/Manual/packages.html#PkgSpatialSearchingD)

-   The kd-tree can now be built in parallel: [`CGAL::Kd_tree::build()`](https://doc.cgal.org/5.1/Spatial_searching/classCGAL_1_1Kd__tree.html#a8559dbe4d7136fbc8ebab5ee290cbe06)
    is given an optional template parameter `ConcurrencyTag` (default
    value remains [`CGAL::Sequential_tag`](https://doc.cgal.org/5.1/STL_Extension/structCGAL_1_1Sequential__tag.html)
    for backward compatibility).
-   Improved the performance of the kd-tree in some cases:
    - Not storing the points coordinates inside the tree usually
      generates a lot of cache misses, leading to non-optimal
      performance. This is the case for example
      when indices are stored inside the tree, or to a lesser extent when the points
      coordinates are stored in a dynamically allocated array (e.g., [`Epick_d`](https://doc.cgal.org/5.1/Kernel_d/structCGAL_1_1Epick__d.html)
      with dynamic dimension) &mdash; we says "to a lesser extent" because the points
      are re-created by the kd-tree in a cache-friendly order after its construction,
      so the coordinates are more likely to be stored in a near-optimal order
      on the heap.
      In these cases, the new `EnablePointsCache` template parameter of the
      [`CGAL::Kd_tree`](https://doc.cgal.org/5.1/Spatial_searching/classCGAL_1_1Kd__tree.html)
      class can be set to `CGAL::Tag_true`. The points coordinates
      will then be cached in an optimal way. This will increase memory
      consumption but provides better search performance. See the updated
      [`GeneralDistance`](https://doc.cgal.org/5.1/Spatial_searching/classGeneralDistance.html)
      and [`FuzzyQueryItem`](https://doc.cgal.org/5.1/Spatial_searching/classFuzzyQueryItem.html)
      concepts for additional requirements when using such a cache.
    - In most cases (e.g., Euclidean distance), the distance computation
      algorithm knows before its end that the distance will be greater
      than or equal to some given value. This is used in the (orthogonal)
      k-NN search to interrupt some distance computations before its end,
      saving precious milliseconds, in particular in medium-to-high dimension.

### [Intersecting Sequences of dD Iso-oriented Boxes](https://doc.cgal.org/5.1/Manual/packages.html#PkgBoxIntersectionD)

-   Added parallel versions of the functions
    [`CGAL::box_intersection_d()`](https://doc.cgal.org/5.1/Box_intersection_d/group__PkgBoxIntersectionD__box__intersection__d.html)
    and [`CGAL::box_self_intersection_d()`](https://doc.cgal.org/5.1/Box_intersection_d/group__PkgBoxIntersectionD__box__self__intersection__d.html).

### [Spatial Sorting](https://doc.cgal.org/5.1/Manual/packages.html#PkgSpatialSorting)

-   Added parallel versions of the functions
    [`CGAL::hilbert_sort()`](https://doc.cgal.org/5.1/Spatial_sorting/group__PkgSpatialSortingFunctions.html#ga9da67204747ac19dff65f9c9ff2fca9e)
    and [`CGAL::spatial_sort()`](https://doc.cgal.org/5.1/Spatial_sorting/group__PkgSpatialSortingFunctions.html#ga7c597c11a3b3859234ff68526cead84d)
    in 2D and 3D when the median policy is used.
    The parallel versions use up to four threads in 2D, and up to eight threads in 3D.

### [3D Convex Hulls](https://doc.cgal.org/5.1/Manual/packages.html#PkgConvexHull3)

-   A new overload for [`CGAL::convex_hull_3()`](https://doc.cgal.org/5.1/Convex_hull_3/group__PkgConvexHull3Functions.html#gaa02a3013808fc9a2e5e2f42b9fde8e30)
    that takes a model of [`VertexListGraph`](https://doc.cgal.org/5.1/BGL/classVertexListGraph.html) has been added.
-   The long-deprecated function `CGAL::convex_hull_3_to_polyhedron_3()` has been removed.
    The function [`CGAL::convex_hull_3_to_face_graph()`](https://doc.cgal.org/5.1/Convex_hull_3/group__PkgConvexHull3Functions.html#ga2750f7f197588ed643679835c748c671)
    should be used instead.

### [Polygon Mesh Processing](https://doc.cgal.org/5.1/Manual/packages.html#PkgPolygonMeshProcessing)

-   Added the function [`CGAL::Polygon_mesh_processing::volume_connected_component()`](https://doc.cgal.org/5.1/Polygon_mesh_processing/group__PMP__orientation__grp.html#ga133e58280959c152770525f27bb42b91),
    which can be used to get information about the nesting of the connected components of a given triangle mesh and about
    the volumes defined.
-   Added the function [`CGAL::Polygon_mesh_processing::remove_connected_components_of_negligible_size()`](https://doc.cgal.org/5.1/Polygon_mesh_processing/group__PMP__repairing__grp.html#gac544fcaba1d59d330a3a1536caff392a),
    which can be used to remove connected components whose area or volume is under a certain threshold.
    Area and volume thresholds are either specified by the user or deduced from the bounding box of the mesh.
-   Added a new named parameter for [`CGAL::Polygon_mesh_processing::keep_large_connected_components()`](https://doc.cgal.org/5.1/Polygon_mesh_processing/group__keep__connected__components__grp.html#ga48e7b3e6922ee78cf8ce801e3e325d9a)
    and [`CGAL::Polygon_mesh_processing::remove_connected_components_of_negligible_size`](https://doc.cgal.org/5.1/Polygon_mesh_processing/group__PMP__repairing__grp.html#gac544fcaba1d59d330a3a1536caff392a),
    which can be used to perform a dry run of the operation, meaning that the function will return the number of connected
    components that would be removed with the specified threshold, but without actually removing them.
-   Added the function [`CGAL::Polygon_mesh_processing::split()`](https://doc.cgal.org/5.1/Polygon_mesh_processing/group__PMP__corefinement__grp.html#gaa491feee9e41f725332bea0ea1215578),
    which can be used to split meshes along a mesh or a plane.
-   Added the function [`CGAL::Polygon_mesh_processing::split_connected_components()`](https://doc.cgal.org/5.1/Polygon_mesh_processing/group__keep__connected__components__grp.html#ga9ddd1e4b915a4232b1ce5611985302aa)
    to split a single mesh containing several connected components into several meshes containing one connected component.
-   Added the functions [`CGAL::Polygon_mesh_processing::merge_reversible_connected_components()`](https://doc.cgal.org/5.1/Polygon_mesh_processing/group__PMP__orientation__grp.html#gae25c1198a89c53d5df2f29dd57fda5ca),
    [`CGAL::Polygon_mesh_processing::duplicate_non_manifold_edges_in_polygon_soup()`](https://doc.cgal.org/5.1/Polygon_mesh_processing/group__PMP__orientation__grp.html#ga2aa4f7b500dc51d1fc4747705a050946),
    and [`CGAL::Polygon_mesh_processing::orient_triangle_soup_with_reference_triangle_mesh()`](https://doc.cgal.org/5.1/Polygon_mesh_processing/group__PMP__orientation__grp.html#ga31779672b3afd660664fc9a6c4fdf74d),
    which can be helpful when repairing a polygon soup.
-   Added the function [`CGAL::Polygon_mesh_processing::sample_triangle_soup()`](https://doc.cgal.org/5.1/Polygon_mesh_processing/group__PMP__distance__grp.html#gac7af41d13bf1a7c30852be266ac81db5),
    which generates points on a triangle soup surface.
-   Added parallel versions of the functions [`CGAL::Polygon_mesh_processing::does_self_intersect()`](https://doc.cgal.org/5.1/Polygon_mesh_processing/group__PMP__intersection__grp.html#gad9fe5d8b433545b69154f43935a11a3b)
    and [`CGAL::Polygon_mesh_processing::self_intersections()`](https://doc.cgal.org/5.1/Polygon_mesh_processing/group__PMP__intersection__grp.html#gaf19c80ec12cbff7ebe9e69453f1d40b8).
-   The function [`CGAL::Polygon_mesh_processing::stitch_borders()`](https://doc.cgal.org/5.1/Polygon_mesh_processing/group__PMP__repairing__grp.html#ga8ae4352e67d2b099994ac8990c13bd41)
    now returns the number of halfedge pairs that were stitched.
-   Added the function [`CGAL::Polygon_mesh_processing::polygon_mesh_to_polygon_soup()`](https://doc.cgal.org/5.1/Polygon_mesh_processing/group__PMP__repairing__grp.html#ga76648a509409ff3c3ad3f71eff8ce9d9).
-   The function [`CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh`](https://doc.cgal.org/5.1/Polygon_mesh_processing/group__PMP__repairing__grp.html#ga0dec58e8a0112791f72ebbe77bac074b)
    now allows passing a point map (for the point range) and a vertex point map (for the polygon mesh) via named parameters.

### [Point Set Processing](https://doc.cgal.org/5.1/Manual/packages.html#PkgPointSetProcessing3)

-   **Breaking change:** [`CGAL::remove_outliers()`](https://doc.cgal.org/5.1/Point_set_processing_3/group__PkgPointSetProcessing3Algorithms.html#ga1ab1dcee59caadde50572c5a504cc41a)
    has been parallelized and thus has a new template parameter `ConcurrencyTag`.
    To update your code simply add as first template parameter `CGAL::Sequential_tag` or `CGAL::Parallel_tag`
    when calling this function.
-   Add a function [`CGAL::cluster_point_set()`](https://doc.cgal.org/5.1/Point_set_processing_3/group__PkgPointSetProcessing3Algorithms.html#gafee41d60b5a257ae034e9157d0af8e46)
    that segments a point cloud into connected components based on a distance threshold.
-   Added wrapper functions for registration:
    - [`CGAL::OpenGR::compute_registration_transformation()`](https://doc.cgal.org/5.1/Point_set_processing_3/group__PkgPointSetProcessing3Algorithms.html#gab81663c718960780ddb176aad845e8cd),
      which computes the registration transformation for two point sets using the Super4PCS algorithm
      implemented in the third party library [OpenGR](https://storm-irit.github.io/OpenGR/index.html).
    - [`CGAL::OpenGR::register_point_sets()`](https://doc.cgal.org/5.1/Point_set_processing_3/group__PkgPointSetProcessing3Algorithms.html#ga6194087f512e4e23dd945a9364d0931d),
      which computes the registration transformation for two point sets using the Super4PCS algorithm
      implemented in the third party library [OpenGR](https://storm-irit.github.io/OpenGR/index.html),
      and registers the points sets by transforming the data point set using the computed transformation.
    - [`CGAL::pointmatcher::compute_registration_transformation()`](https://doc.cgal.org/5.1/Point_set_processing_3/group__PkgPointSetProcessing3Algorithms.html#gaf75af5c1634fa83fa05a33e95570b127)
      computes the registration transformation for two point sets using ICP algorithm implemented
      in the third party library [libpointmatcher](https://github.com/ethz-asl/libpointmatcher).
    - [`CGAL::pointmatcher::register_point_sets()`](https://doc.cgal.org/5.1/Point_set_processing_3/group__PkgPointSetProcessing3Algorithms.html#gaa222278e20a3ce41930d37326cd54ef9),
      which computes the registration transformation for two point sets using ICP algorithm implemented
      in the third party library [libpointmatcher](https://github.com/ethz-asl/libpointmatcher), and registers
      the points sets by transforming the data point set using the computed transformation.

### [2D Triangulations](https://doc.cgal.org/5.1/Manual/packages.html#PkgTriangulation2)

-   To fix an inconsistency between code and documentation and to clarify which types of intersections
    are truly allowed in constrained Delaunay triangulations, the tag [`CGAL::No_intersection_tag`](https://doc.cgal.org/5.1/Triangulation_2/structCGAL_1_1No__intersection__tag.html)
    has been deprecated in favor of two new tags: [`CGAL::No_constraint_intersection_tag`](https://doc.cgal.org/5.1/Triangulation_2/structCGAL_1_1No__constraint__intersection__tag.html)
    and [`CGAL::No_constraint_intersection_requiring_constructions_tag`](https://doc.cgal.org/5.1/Triangulation_2/structCGAL_1_1No__constraint__intersection__requiring__constructions__tag.html).
    The latter is equivalent to the now-deprecated `CGAL::No_intersection_tag`, and allows constraints
    to intersect as long as no new point has to be created to represent that intersection (for example,
    the intersection of two constraint segments in a 'T'-like junction is an existing point
    and as such does not require any new construction). The former tag, `CGAL::No_constraint_intersection_tag`,
    does not allow any intersection, except for the configuration of two constraints having a single
    common endpoints, for convenience.
-   Added the function [`CGAL::split_subconstraint_graph_into_constraints()`](https://doc.cgal.org/5.1/Triangulation_2/classCGAL_1_1Constrained__triangulation__plus__2.html#adea77f5db5cd4dfae302e4502f1caa85)
    to [`Constrained_triangulation_plus_2`](https://doc.cgal.org/5.1/Triangulation_2/classCGAL_1_1Constrained__triangulation__plus__2.html) to initialize the constraints
    from a soup of disconnected segments that should first be split into polylines.

### [3D Triangulations](https://doc.cgal.org/5.1/Manual/packages.html#PkgTriangulation3)

-   The member function [`CGAL::Triangulation_3::file_input()`](https://doc.cgal.org/5.1/Triangulation_3/group__PkgIOTriangulation3.html#gadd94d0613e2dd9cdd2e88d2c74d5b1c8)
    have been added. It allows to load a [`CGAL::Triangulation_3`](https://doc.cgal.org/5.1/Triangulation_3/classCGAL_1_1Triangulation__3.html)
    from an input stream, using functors to create vertices and cells.

### [3D Triangulation Data Structure](https://doc.cgal.org/5.1/Manual/packages.html#PkgTDS3)

-   The member function [`CGAL::TDS_3::file_input()`](https://doc.cgal.org/5.1/TDS_3/group__PkgIOTDS3.html#ga381446a02a9240cc83e79c48b37cd119)
    have been added. It allows to load a [`CGAL::Triangulation_data_structure_3`](https://doc.cgal.org/5.1/TDS_3/classCGAL_1_1Triangulation__data__structure__3.html)
    from an input stream, using functors to create vertices and cells.

### [Surface Mesh Simplification](https://doc.cgal.org/5.1/Manual/packages.html#PkgSurfaceMeshSimplification)

-   Added a [new simplification method](https://doc.cgal.org/5.1/Surface_mesh_simplification/classCGAL_1_1Surface__mesh__simplification_1_1GarlandHeckbert__policies.html)
    based on the quadric error defined by Garland and Heckbert.
-   The concept `EdgeProfile` has been removed. This concept was not actually in use as the CGAL-provided model [`CGAL::Edge_profile`](https://doc.cgal.org/5.1/Surface_mesh_simplification/classCGAL_1_1Surface__mesh__simplification_1_1Edge__profile.html)
    was imposed to the user. Other concepts have been clarified to reflect the fact that the API uses this particular class.

### [STL Extensions for CGAL](https://doc.cgal.org/5.1/Manual/packages.html#PkgSTLExtension)

-   Added a new concurrency tag: [`CGAL::Parallel_if_available_tag`](https://doc.cgal.org/5.1/STL_Extension/structCGAL_1_1Parallel__if__available__tag.html).
    This tag is a convenience typedef to [`CGAL::Parallel_tag`](https://doc.cgal.org/5.1/STL_Extension/structCGAL_1_1Parallel__tag.html)
    if the third party library TBB has been found and linked with, and to
    [`CGAL::Sequential_tag`](https://doc.cgal.org/5.1/STL_Extension/structCGAL_1_1Sequential__tag.html) otherwise.


[Release 5.0](https://github.com/CGAL/cgal/releases/tag/releases%2FCGAL-5.0)
-----------

Release date: November 2019

### General changes

- CGAL 5.0 is the first release of CGAL that requires a C++ compiler
  with the support of C++14 or later. The new list of supported
  compilers is:
  - Visual C++ 14.0 (from Visual Studio 2015 Update 3) or later,
  - Gnu g++ 6.3 or later (on Linux or MacOS),
  - LLVM Clang version 8.0 or later (on Linux or MacOS), and
  - Apple Clang compiler versions 7.0.2 and 10.0.1 (on MacOS).
- Since CGAL 4.9, CGAL can be used as a header-only library, with
  dependencies. Since CGAL 5.0, that is now the default, unless
  specified differently in the (optional) CMake configuration.
- The section "Getting Started with CGAL" of the documentation has
  been updated and reorganized.
- The minimal version of Boost is now 1.57.0.


### [Polygonal Surface Reconstruction](https://doc.cgal.org/5.0/Manual/packages.html#PkgPolygonalSurfaceReconstruction) (new package)

 -   This package provides a method for piecewise planar object reconstruction from point clouds.
     The method takes as input an unordered point set sampled from a piecewise planar object
     and outputs a compact and watertight surface mesh interpolating the input point set.
     The method assumes that all necessary major planes are provided (or can be extracted from
     the input point set using the shape detection method described in Point Set Shape Detection,
     or any other alternative methods).The method can handle arbitrary piecewise planar objects
     and is capable of recovering sharp features and is robust to noise and outliers. See also
     the associated [blog entry](https://www.cgal.org/2019/08/05/Polygonal_surface_reconstruction/).

### [Shape Detection](https://doc.cgal.org/5.0/Manual/packages.html#PkgShapeDetection) (major changes)
 -   **Breaking change:** The concept `ShapeDetectionTraits` has been renamed to [`EfficientRANSACTraits`](https://doc.cgal.org/5.0/Shape_detection/classEfficientRANSACTraits.html).
 -   **Breaking change:** The `Shape_detection_3` namespace has been renamed to [`Shape_detection`](https://doc.cgal.org/5.0/Shape_detection/annotated.html).
 -   Added a new, generic implementation of region growing. This enables for example applying region growing to inputs such as 2D and 3D point sets,
     or models of the [`FaceGraph`](https://doc.cgal.org/5.0/BGL/classFaceGraph.html) concept. Learn more about this new algorithm with this [blog entry](https://www.cgal.org/2019/07/30/Shape_detection/).

### [dD Geometry Kernel](https://doc.cgal.org/5.0/Manual/packages.html#PkgKernelD)
 -   A new exact kernel, [`Epeck_d`](https://doc.cgal.org/5.0/Kernel_d/structCGAL_1_1Epeck__d.html), is now available.

### [2D and 3D Linear Geometry Kernel](https://doc.cgal.org/5.0/Manual/packages.html#PkgKernel23)
 -   Added a new concept, [`ComputeApproximateAngle_3`](https://doc.cgal.org/5.0/Kernel_23/classKernel_1_1ComputeApproximateAngle__3.html),
     to the 3D Kernel concepts to compute the approximate angle between two 3D vectors. Corresponding functors
     in the model ([`Compute_approximate_angle_3`](https://doc.cgal.org/5.0/Kernel_23/classKernel.html#a183c9ac358a4ccddc04e680f8ed16c0b))
     and free function ([`approximate_angle`](https://doc.cgal.org/5.0/Kernel_23/group__approximate__angle__grp.html))
     have also been added.
 -   The following objects are now hashable and thus trivially usable
     with [`std::unordered_set`](https://en.cppreference.com/w/cpp/container/unordered_set)
     and [`std::unordered_map`](https://en.cppreference.com/w/cpp/header/unordered_map):
     `CGAL::Aff_transformation_2`, `CGAL::Aff_transformation_3`,
     `CGAL::Bbox_2`, `CGAL::Bbox_3`, `CGAL::Circle_2`,
     `CGAL::Iso_cuboid_3`, `CGAL::Iso_rectangle_2`, `CGAL::Point_2`,
     `CGAL::Point_3`, `CGAL::Segment_2`, `CGAL::Segment_3`,
     `CGAL::Sphere_3`, `CGAL::Vector_2`, `CGAL::Vector_3`,
     `CGAL::Weighted_point_2` and `CGAL::Weighted_point_3`.

### [Polygon Mesh Processing](https://doc.cgal.org/5.0/Manual/packages.html#PkgPolygonMeshProcessing)
 -   Introduced a [wide range of new functions](https://doc.cgal.org/5.0/Polygon_mesh_processing/index.html#title36)
     related to location of queries on a triangle mesh,
     such as [`CGAL::Polygon_mesh_processing::locate(Point, Mesh)`](https://doc.cgal.org/5.0/Polygon_mesh_processing/group__PMP__locate__grp.html#gada09bd8740ba69ead9deca597d53cf15).
     The location of a point on a triangle mesh is expressed as the pair of a face and the barycentric
     coordinates of the point in this face, enabling robust manipulation of locations
     (for example, intersections of two 3D segments living within the same face).
 -   Added the mesh smoothing function [`smooth_mesh()`](https://doc.cgal.org/5.0/Polygon_mesh_processing/group__PMP__meshing__grp.html#gaa0551d546f6ab2cd9402bea12d8332a3),
     which can be used to improve the quality of triangle elements based on various geometric characteristics.
 -   Added the shape smoothing function [`smooth_shape()`](https://doc.cgal.org/5.0/Polygon_mesh_processing/group__PMP__meshing__grp.html#gaaa083ec78bcecf351e04d1bbf460b4a2),
     which can be used to smooth the surface of a triangle mesh, using the mean curvature flow to perform noise removal.
     (See also the new entry in the [User Manual](https://doc.cgal.org/5.0/Polygon_mesh_processing/index.html#title8))
 -   Added the function [`CGAL::Polygon_mesh_processing::centroid()`](https://doc.cgal.org/5.0/Polygon_mesh_processing/group__measure__grp.html#ga6da5119ce2c50729fda11a90ae7fb9ba),
     which computes the centroid of a closed triangle mesh.
 -   Added the functions [`CGAL::Polygon_mesh_processing::stitch_boundary_cycle()`](https://doc.cgal.org/5.0/Polygon_mesh_processing/group__PMP__repairing__grp.html#ga9c12c4878c08a117b3733bb45f1a34cf)
     and [`CGAL::Polygon_mesh_processing::stitch_boundary_cycles()`](https://doc.cgal.org/5.0/Polygon_mesh_processing/group__PMP__repairing__grp.html#ga24d5ae37f62064b3fc576ba48a4ccc63),
     which can be used to try and merge together geometrically compatible but combinatorially different halfedges
     that belong to the same boundary cycle.
 -   It is now possible to pass a face-size property map to [`CGAL::Polygon_mesh_processing::keep_large_connected_components()`](https://doc.cgal.org/5.0/Polygon_mesh_processing/group__keep__connected__components__grp.html#ga48e7b3e6922ee78cf8ce801e3e325d9a)
     and [`CGAL::Polygon_mesh_processing::keep_largest_connected_components()`](https://doc.cgal.org/5.0/Polygon_mesh_processing/group__keep__connected__components__grp.html#ga68c6c29dfc6a26a6a2f8befe6944f19d), enabling users to define
     how the size of a face is computed (the size of the connected component is the sum of the sizes of its faces).
     If no property map is passed, the behavior is unchanged to previous versions: the size
     of a connected component is the number of faces it contains.
 -   Added the function [`CGAL::Polygon_mesh_processing::non_manifold_vertices()`](https://doc.cgal.org/5.0/Polygon_mesh_processing/group__PMP__repairing__grp.html#ga36098d2415efd0604b7b996163bc22db),
     which can be used to collect all the non-manifold vertices (i.e. pinched vertices,
     or vertices appearing in multiple umbrellas) of a mesh.

### [3D Point Set](https://doc.cgal.org/5.0/Manual/packages.html#PkgPointSet3)
 -   The [PLY IO functions](https://doc.cgal.org/5.0/Point_set_3/group__PkgPointSet3IO.html) now take an additional optional parameter to
     read/write comments from/in the PLY header.

### [Point Set Processing](https://doc.cgal.org/5.0/Manual/packages.html#PkgPointSetProcessing3)
 -   **Breaking change**: the API using iterators and overloads for optional parameters (deprecated since
     CGAL 4.12) has been removed. The current (and now only) API uses ranges and Named Parameters.
 -   Added the possibility to use the named parameter
     [`neighbor_radius`](https://doc.cgal.org/5.0/Point_set_processing_3/group__psp__namedparameters.html#PSP_neighbor_radius)
     to use spherical neighbor queries instead of K-nearest neighbors queries for the following functions:
     [`CGAL::bilateral_smooth_point_set()`](https://doc.cgal.org/5.0/Point_set_processing_3/group__PkgPointSetProcessing3Algorithms.html#ga4f82723e2f0bb33f3677e29e0208a256),
     [`CGAL::jet_estimate_normals()`](https://doc.cgal.org/5.0/Point_set_processing_3/group__PkgPointSetProcessing3Algorithms.html#ga0cd0f87de690d4edf82740e856efa491),
     [`CGAL::jet_smooth_point_set()`](https://doc.cgal.org/5.0/Point_set_processing_3/group__PkgPointSetProcessing3Algorithms.html#ga549402c0a8a8b6b71875181e93961521),
     [`CGAL::mst_orient_normals()`](https://doc.cgal.org/5.0/Point_set_processing_3/group__PkgPointSetProcessing3Algorithms.html#ga50c98d5c5ae5535bce6f32eddbd03f33),
     [`CGAL::pca_estimate_normals()`](https://doc.cgal.org/5.0/Point_set_processing_3/group__PkgPointSetProcessing3Algorithms.html#ga8c642da96a025ab32445aeb6cc219b0b) and
     [`CGAL::remove_outliers()`](https://doc.cgal.org/5.0/Point_set_processing_3/group__PkgPointSetProcessing3Algorithms.html#gafd0b5a21ec5042e4bca09cb43f1847f9).

### [2D Triangulations](https://doc.cgal.org/5.0/Manual/packages.html#PkgTriangulation2)
 -   **Breaking change**: Removed the deprecated functions `CGAL::Constrained_triangulation_plus_2::
     vertices_in_constraint_{begin/end}(Vertex_handle va, Vertex_handle vb) const;`,
     and `CGAL::Constrained_triangulation_plus_2::remove_constraint(Vertex_handle va, Vertex_handle vb)`,
     that is a pair of vertex handles is no longer a key for a polyline constraint.
     Users must use a version prior to 5.0 if they need this functionality.
 -   **Breaking change**: Removed the deprecated classes `CGAL::Regular_triangulation_euclidean_traits_2`,
     `CGAL::Regular_triangulation_filtered_traits_2`. Users must use a version prior to 5.0 if they need these classes.
 -   **Breaking change**: The [graph traits](https://doc.cgal.org/5.0/BGL/group__PkgBGLTraits.html) enabling CGAL's 2D triangulations to be used as a parameter
     for any graph-based algorithm of CGAL (or boost) have been improved to fully model the [`FaceGraph`](https://doc.cgal.org/5.0/BGL/classFaceGraph.html) concept.
     In addition, only the finite simplicies (those not incident to the infinite vertex) of the 2D triangulations
     are now visible through this scope. The complete triangulation can still be accessed as a graph,
     by using the graph traits of the underlying triangulation data structure (usually,
     [`CGAL::Triangulation_data_structure_2`](https://doc.cgal.org/5.0/TDS_2/classCGAL_1_1Triangulation__data__structure__2.html)).
 -   **Breaking change**: The `insert()` function
     of
     [`CGAL::Triangulation_2`](https://doc.cgal.org/5.0/Triangulation_2/classCGAL_1_1Triangulation__2.html)
     which takes a range of points as argument is now guaranteed to
     insert the points following the order of `InputIterator`.  Note
     that this change only affects the base class `Triangulation_2`
     and not any derived class, such as `Delaunay_triangulation_2`.
-   Added a new [constructor](https://doc.cgal.org/5.0/Triangulation_2/classCGAL_1_1Triangulation__2.html#a6cfa7d3aaa375a25d217858b49e2eb07=)
     and [`insert()`](https://doc.cgal.org/5.0/Triangulation_2/classCGAL_1_1Triangulation__2.html#ac5e9bc8adef80dc01a0b31c2d0234545)
     function to [`CGAL::Triangulation_2`](https://doc.cgal.org/5.0/Triangulation_2/classCGAL_1_1Triangulation__2.html)
     that takes a range of points with info.
 -   Introduced a new face base class, [`Triangulation_face_base_with_id_2`](https://doc.cgal.org/5.0/BGL/classCGAL_1_1Triangulation__face__base__with__id__2.html)
     which enables storing user-defined integer IDs in the face of any 2D triangulation, a precondition to use some
     BGL algorithms.
 -   Added range types and functions that return ranges, for example for all vertices, enabling the use of `C++11` `for`-loops.
     See [this new example](https://doc.cgal.org/5.0/Triangulation_2/Triangulation_2_2for_loop_2_8cpp-example.html) for a usage demonstration.

### [3D Triangulations](https://doc.cgal.org/5.0/Manual/packages.html#PkgTriangulation3)
 -   **Breaking change**: The [constructor](https://doc.cgal.org/5.0/Triangulation_3/classCGAL_1_1Triangulation__3.html#a63f67cf6aaadcee14318cf56a36d247a)
     and the [`insert()`](https://doc.cgal.org/5.0/Triangulation_3/classCGAL_1_1Triangulation__3.html#ad3353128386bbb51f79d0263e7f67337)
     function of [`CGAL::Triangulation_3`](https://doc.cgal.org/5.0/Triangulation_3/classCGAL_1_1Triangulation__3.html)
     which take a range of points as argument are now guaranteed to
     insert the points following the order of `InputIterator`. Note
     that this change only affects the base class `Triangulation_3`
     and not any derived class, such as `Delaunay_triangulation_3`.
 -   Added constructor and [`insert()`](https://doc.cgal.org/5.0/Triangulation_3/classCGAL_1_1Triangulation__3.html#a8aa85f88733d30aa3ec5385538e13ace)
     function to `CGAL::Triangulation_3` that takes a range of points with info.
 -   Added range types and functions that return ranges, for example for all vertices, which enables to use C++11 for-loops.
     See [this new example](https://doc.cgal.org/5.0/Triangulation_3/Triangulation_3_2for_loop_8cpp-example.html) for a usage demonstration.

### [Surface Mesh](https://doc.cgal.org/5.0/Manual/packages.html#PkgSurfaceMesh)
 -   Introduced new functions to read and write using the PLY format,
     [`CGAL::read_ply()`](https://doc.cgal.org/5.0/Surface_mesh/group__PkgSurface__mesh.html#ga42f6ad486ddab74e13d3dc53f511c343)
     and [`CGAL::write_ply()`](https://doc.cgal.org/5.0/Surface_mesh/group__PkgSurface__mesh.html#ga77bbb79d449c981895eedb6c3c23bd14),
     enabling users to save and load additional property maps of the surface mesh.

###  [CGAL and Solvers](https://doc.cgal.org/5.0/Manual/packages.html#PkgSolverInterface)
 -   Added [concepts](https://doc.cgal.org/5.0/Solver_interface/group__PkgSolverInterfaceConcepts.html)
     and [models](https://doc.cgal.org/5.0/Solver_interface/group__PkgSolverInterfaceRef.html)
     for solving Mixed Integer Programming (MIP) problems with or without constraints.

### [3D Boolean Operations on Nef Polyhedra](https://doc.cgal.org/5.0/Manual/packages.html#PkgNef3)
 -   Added a function to convert a Nef_polyhedron_3 to a polygon soup: [`CGAL::convert_nef_polyhedron_to_polygon_soup()`](https://doc.cgal.org/5.0/Nef_3/group__PkgNef3IOFunctions.html#ga28a9eb4da0cd6153f0c16f7f9eaf6665)

### [IO Streams](https://doc.cgal.org/5.0/Manual/packages.html#PkgStreamSupport)
- **Breaking change:** The API of [`CGAL::Color`](https://doc.cgal.org/5.0/Stream_support/classCGAL_1_1Color.html) has been cleaned up.
- Added new functions to support some parts of the WKT file format:
    * [`CGAL::read_WKT()`](https://doc.cgal.org/5.0/Stream_support/group__PkgStreamSupportRef.html#gad2872abfe6fcf17d705d38567fdd6248)
    * [`CGAL::read_point_WKT()`](https://doc.cgal.org/5.0/Stream_support/group__PkgStreamSupportRef.html#gadbd2705b183e467507abd2f167446eba)
    * [`CGAL::read_multi_point_WKT()`](https://doc.cgal.org/5.0/Stream_support/group__PkgStreamSupportRef.html#ga4fb72e49a1fd385bbed35ea20297aa8d)
    * [`CGAL::read_linestring_WKT()`](https://doc.cgal.org/5.0/Stream_support/group__PkgStreamSupportRef.html#gaaa236308b9da5dbf217ef281fdb55de4)
    * [`CGAL::read_multi_linestring_WKT()`](https://doc.cgal.org/5.0/Stream_support/group__PkgStreamSupportRef.html#gad6046c7f9d36512b8a014be82c1e2220)
    * [`CGAL::read_polygon_WKT()`](https://doc.cgal.org/5.0/Stream_support/group__PkgStreamSupportRef.html#gaa36ccd3ac4b3fe3e3fd8a76715c56b9a)
    * [`CGAL::read_multi_polygon_WKT()`](https://doc.cgal.org/5.0/Stream_support/group__PkgStreamSupportRef.html#ga4ceaa71b9cb3b3f7984bed19afff6fc6)
    * [`CGAL::write_point_WKT()`](https://doc.cgal.org/5.0/Stream_support/group__PkgStreamSupportRef.html#gab1a2d277b43c218bf128a2056eb53ced)
    * [`CGAL::write_polygon_WKT()`](https://doc.cgal.org/5.0/Stream_support/group__PkgStreamSupportRef.html#gab5365a4726893aa4f51739ede63f5a09)
    * [`CGAL::write_linestring_WKT()`](https://doc.cgal.org/5.0/Stream_support/group__PkgStreamSupportRef.html#gaa37ed77d1a01567b93c872a48198efa6)
    * [`CGAL::write_multi_point_WKT()`](https://doc.cgal.org/5.0/Stream_support/group__PkgStreamSupportRef.html#ga98de4b4e5cccb370febe5daf66bb582d)
    * [`CGAL::write_multi_polygon_WKT()`](https://doc.cgal.org/5.0/Stream_support/group__PkgStreamSupportRef.html#ga4ded40ab50f57e0b410640e28964935e)
    * [`CGAL::write_multi_linestring_WKT()`](https://doc.cgal.org/5.0/Stream_support/group__PkgStreamSupportRef.html#ga219987f7a9c0b871c1733aa0c38f26b3)


Release 4.14
------------

Release date: March 2019

### 2D Periodic Hyperbolic Triangulations (new package)

 -   This package allows the computation of Delaunay triangulations of
     the Bolza surface.  The Bolza surface is the most symmetric
     hyperbolic surface of genus 2. Its fundamental domain is the
     regular hyperbolic octagon with angles π/4 centered at the origin
     of the Poincaré disk. Triangulations of the Bolza surface can be
     seen as triangulations of the hyperbolic plane that are periodic
     in the four directions defined by the sides of this regular
     octagon.

### 2D Hyperbolic Triangulations (new package)

 -   This package allows the computation of Delaunay Triangulations of
     sets of points in the Poincaré disk, which is one of the
     conformal models for the hyperbolic plane.

### The Heat Method (new package)

-   This package provides an algorithm that solves the single- or
    multiple-source shortest path problem by returning an
    approximation of the geodesic distance for all vertices of a
    triangle mesh to the closest vertex in a given set of source
    vertices.

### Triangulated Surface Mesh Approximation (new package)

-   This package implements the Variational Shape Approximation method
    to approximate an input surface triangle mesh by a simpler surface
    triangle mesh.

### Polygon Mesh Processing package

-   Added the following new functions to detect and repair issues in
    polygon soups:
    - `CGAL::Polygon_mesh_processing::remove_isolated_points_in_polygon_soup()`,
       which detects and removes points that are not used in any
       polygon of the soup.
    - `CGAL::Polygon_mesh_processing::merge_duplicate_points_in_polygon_soup()`,
       which detects and merges points that share the same geometric
       position.
    - `CGAL::Polygon_mesh_processing::merge_duplicate_polygons_in_polygon_soup()`,
       which detects and merges polygons that are identical.
    - `CGAL::Polygon_mesh_processing::repair_polygon_soup()`, which
       applies a number of repairing steps (a subset of which are the
       functions above) to clean and repair a polygon soup.

-   Added the following new functions to detect and repair
    degeneracies in polygon meshes:
    - `CGAL::Polygon_mesh_processing::degenerate_edges()`
    - `CGAL::Polygon_mesh_processing::degenerate_faces()`
    - `CGAL::Polygon_mesh_processing::is_non_manifold_vertex()`
    - `CGAL::Polygon_mesh_processing::is_degenerate_triangle_face()`
    - `CGAL::Polygon_mesh_processing::is_degenerate_edge()`
    - `CGAL::Polygon_mesh_processing::is_needle_triangle_face()`
    - `CGAL::Polygon_mesh_processing::is_cap_triangle_face()`
    - `CGAL::Polygon_mesh_processing::duplicate_non_manifold_vertices()`
    - `CGAL::Polygon_mesh_processing::extract_boundary_cycles()`
    - `CGAL::Polygon_mesh_processing::merge_duplicated_vertices_in_boundary_cycle()`
    - `CGAL::Polygon_mesh_processing::merge_duplicated_vertices_in_boundary_cycles()`

-   Added the class `CGAL::Rigid_triangle_mesh_collision_detection` to
    detect intersections between meshes and volumes undergoing affine
    transformations.

### Regularized Boolean Set Operations in 2D package

-   Fixed the validation of orientation of relative simple polygons.

### Point Set Processing

-   `CGAL::mst_orient_normals()` can now be called with a set of
    user-selected seed points that are known to be already oriented. A
    new optional named parameter `point_is_constrained_map` is added
    for this purpose. The original behavior (using one unique and
    automatically selected seed) is kept if this parameter is not
    used.

### Classification

-   Added a new experimental classifier
    `TensorFlow::Neural_network_classifier`.

-   For uniformity, `ETHZ_random_forest_classifier` is renamed
    `ETHZ::Random_forest_classifier` and
    `OpenCV_random_forest_classifier` is renamed
    `OpenCV::Random_forest_classifier`.

-   The training algorithm of `ETHZ::Random_forest_classifier` was
    parallelized.

-   Added a constructor to copy a `ETHZ::Random_forest_classifier`
    using a different data set as input.

-   Added 3 new geometric features, `Height_above`, `Height_below` and
    `Vertical_range`.

### 3D Fast Intersection and Distance Computation

-   The primitives `AABB_face_graph_triangle_primitive` and
    `AABB_halfedge_graph_segment_primitive` now use as `Id` a pair of
    descriptor and graph pointer in the case they are configured to
    deal with a possible different graph per primitive (configuration
    set using a template tag).

### 2D Arrangements

-   Fixed a bug in the surface-sweep framework (`Surface_sweep_2`)
    that ensures that an event is never left without (left or right)
    curves.

-   Fixed a constructor of `Arr_counting_traits.h`. (In particular,
    added missing const of a parameter).

-   Fixed zone computation of a curve in cases where the lexicographic
    smallest end of the curve lies on the parameter space.

-   Implemented missing function object `Compare_x_near_boundary` of
    `Arr_polyline_traits_2`, `Arr_polycurve_traits_2`, and
    `Arr_polycurve_basic_traits_2`.

### 2D and 3D Mesh Generation

-   Added two functions for writing in XML VTK formats:
    - `CGAL::write_vtu()`, that writes a 2D mesh in a `.vtu` file,
    - `CGAL::output_to_vtu()`, that writes a 3D mesh in a `.vtu` file.

### 2D Minkowski Sums

-   Fixed a bug in the function that computed the Minkowski sum using
    the reduced-convolution method. In particular, correctly handled
    the case where one of the summands does not have an outer
    boundary.

### 3D Point Set

-   Added a method `copy_properties()` that allows to copy the
    properties from a point set to another one (without copying the
    content);

-   Added a method `insert(const Point_set&, const Index&)` to copy a
    point along with all its associated properties from another point
    set;

-   `remove()` methods now only invalidate the `end()` iterator
    instead of invalidating all iterators;

-   Added a method `is_removed()` that takes an index as argument;

-   Added a method `cancel_removals()` to restore removed points (if
    no point was inserted since then an garbage was not collected);

-   **Breaking change:** unified API of method `add_normal_map()` with
    `add_property_map()`: it now returns a pair of property map + bool
    (that tells if the property was added) instead of just the
    property map;

-   Added a method `properties_and_types()` in addition to
    `properties()`: this new one returns pairs of `std::string` +
    `std::type_info` in order to also know the type of each property.

### CGAL and the Boost Graph Library (BGL)

-   Added function `write_wrl()` for writing into VRML 2.0 format.
-   Added functions `CGAL::write_vtp()` for writing a triangulated
      face graph in a `.vtp` file (XML VTK format).


Release 4.13
------------

Release date: October 2018

### 3D Periodic Mesh Generation (new package)

-   This package generates 3-dimensional periodic meshes. It computes
    isotropic simplicial meshes for domains described through implicit
    functional boundaries over the flat torus (which can also seen in
    the Euclidean space as a periodic cube). The output is a periodic
    3D mesh of the domain volume and conformal surface meshes for all
    the boundary and subdividing surfaces.  The package is closely
    related to the 3D Mesh Generation package, with similar concepts,
    classes, and API.

### Installation

-   The library `CGAL_Qt5` now contains a fork of the version 2.7.0 of
    `libQGLViewer`.  The corresponding code is in the package
    `GraphicsView`.  The dependency for the external library
    `libQGLViewer` is therefore dropped for all demos.

### General

 -  A new function `CGAL::draw()` is added in the packages Polyhedral
    Surface, Surface Mesh, Linear Cell Complex, 2D Triangulations, and
    3D Triangulations, enabling to draw the corresponding data
    structures.

### 2D and 3D Linear Geometry Kernel

-   An `operator()` that takes a `Ray_3` has been added to the concept
    `ConstructProjectedPoint_3`.

### Convex hull 3

-   Added the function `extreme_points_3()` computing the points on
    the convex hull without underlying connectivity.

-   Added a traits adapter called `Extreme_points_traits_adapter_3`
    that enables the use of the function `extreme_points_3()` on a
    range of keys, each key being associated to 3D point using a
    property map.  This can be used to get the vertices of a mesh that
    are on it convex hull, or the indices of points in a range that
    are on it convex hull.

-   Fix a bug in the computation of the 3D convex hull that was
    leaving extra points within subset of coplanar points that do not
    belong to the minimal convex hull.


### 2D and 3D Triangulations

-   Added a new type of intersection to handle the insertion of
    intersecting constraints in a `Constrained_triangulation_2`.

-   **Breaking change:** The long-deprecated class
    `Triangulation_cell_base_with_circumcenter_3` and its associated
    concept have been removed. Users should use the classes
    `Delaunay_cell_base_with_circumcenter_3` or
    `Regular_cell_base_with_circumcenter_3`, depending on which type
    of triangulation they are using.

-   **Breaking change:** The deprecated functions `mirror_index` and
    `mirror_vertex` of the class `Triangulation_face_base_2` have been
    removed. Users should use the equivalent functions from the class
    `Triangulation_2`.

### 3D Mesh Generation

-   **Breaking change:** The template parameters of the class template
    `Labeled_mesh_domain_3` have been simplified. The three
    constructors of that class template have been replaced by a new
    unique constructor using Boost named parameters. Three new static
    template member functions that act as named constructors have been
    added:
      - `create_gray_image_mesh_domain()`, to create a domain from a 3D
        gray image,
      - `create_labeled_image_mesh_domain()`, to create a domain from a 3D
        labeled image, and
      - `create_implicit_mesh_domain()`, to create a domain from an
        implicit function.

-   The class templates `Implicit_mesh_domain_3`,
    `Gray_image_mesh_domain_3`, and `Labeled_image_mesh_domain_3` are
    now deprecated.

-   **Breaking change:** The headers
    `<CGAL/Mesh_3/Implicit_to_labeled_function_wrapper.h>` and
    `<CGAL/Mesh_3/Labeled_mesh_domain_3.h>`, that were deprecated
    since CGAL 4.5, have been removed.

-   **Breaking change:** the concepts `MeshCellCriteria_3` and
    `MeshFacetCriteria_3` now require the triangulation to be passed
    in their `operator()`.  Models of these concepts that are provided
    by CGAL have been modified accordingly.

-   **Breaking change:** It is no longer possible to use the
    deprecated, pre-CGAL 3.8 specifications in `MeshCellCriteria_3`
    and `MeshFacetCriteria_3` (that is, using `Facet_badness` and
    `Cell_badness` instead of `Is_facet_bad` and `Is_cell_bad`).

-   The concept `MeshFacetCriteria_3` no longer requires the function
    `operator()(Cell_handle c, int i)`.

-   The concept `MeshEdgeCriteria_3` no longer requires the function
    `operator()(const Edge& e)`.

-   The concept `MeshComplexWithFeatures_3InTriangulation_3` no longer
    requires the functions `number_of_edges(Curve_index index)` and
    `number_of_corners(Corner_index index)`.

-   Introduced the concept `MeshTriangulationTraits_3`, which covers
    the needs of the traits class used in `Mesh_3` (and
    `Periodic_3_mesh_3`). The traits class used as template parameter
    of `Mesh_triangulation_3` and `Periodic_3_mesh_triangulation_3`
    must be a model of this concept.

-   Added the function
    `Mesh_domain_with_polyline_features_3::add_corner()`, which allows
    users to add a single corner (that is not incident to any
    polyline) to the mesh complex.

-   **Breaking change**: `CGAL::lloyd_optimize_mesh_3` now depends on
    the _Eigen_ library.

### Polygon Mesh Processing

-   Added a named parameter in stitching functions that allows to
    choose whether the operation should be performed per connected
    component or globally.

-   Added a function, `CGAL::Polygon_mesh_processing::transform()`, to
    apply a transformation to a mesh.

-   Added a named parameter `visitor` in corefinement-related
    functions that makes it possible to pass a visitor to the function
    in order to track the creation of new faces.

-   Added a named parameter `throw_on_self_intersection` in all
    corefinement-related functions, which enables to check for
    self-intersecting faces involved in the intersection before trying
    to corefine the input meshes. This new parameter replaces the
    `bool` parameter in `corefine()`.

-   Added the function `corefine_and_compute_boolean_operations()`,
    which can be used to compute the result of several Boolean
    operations between two volumes at the same time.

-   Added the function `clip()`, which can be used to clip a
    triangulated surface mesh by a plane or a clipping volume.

-   Constrained vertices are now guaranteed to be kept in the mesh
    after calling `isotropic_remeshing()` (and not only the points
    associated to constrained vertices, as it was before).

-   Added a function, `CGAL::Polygon_mesh_processing::extrude_mesh()`,
    to perform an extrusion of an open polygon mesh.

### Estimation of Local Differential Properties of Point-Sampled Surfaces Reference

-   **Breaking change**: `CGAL::Monge_via_jet_fitting` now depends on
    the _Eigen_ library.

### Point Set Processing

-   Added a callback mechanism to the following functions:
    `CGAL::bilateral_smooth_point_set()`,
    `CGAL::compute_average_spacing()`,
    `CGAL::grid_simplify_point_set()`,
    `CGAL::hierarchy_simplify_point_set()`,
    `CGAL::jet_estimate_normals()`, `CGAL::jet_smooth_point_set()`,
    `CGAL::pca_estimate_normals()`, `CGAL::remove_outliers()` and
    `CGAL::wlop_simplify_and_regularize_point_set()`.


### Classification

-   Added data structures to handle classification of Surface Meshes
    and of Clusters.

-   Added public API to compute features in parallel.

-   **Breaking change**: features based on products/divisions of
    eigenvalues are replaced by simple eigenvalue features. Features
    based on statistics on the HSV color channels are replaced by
    simple HSV color channel features.

-   **Breaking change**: the API of
    `CGAL::Classification::Point_set_feature_generator` has been
    simplified.


### Bounding Volumes

-   **Breaking change**: `CGAL::Approximate_min_ellipsoid_d` now
    depends on the _Eigen_ library.

### Interpolation

-   The output of the natural and regular neighbor functions
    (resp. the gradient fitting functions) is no longer restricted to
    a Point/Coordinate pair (resp. Point/Vector pair). Instead, users
    can provide their own functor to format the output as they desire.

-   The interpolation functions can now operate on any combination of
    Type/Coordinate, provided that the values and gradients functors
    can also be evaluated using 'Type'.

    The combination of these two changes allow, for example, to
    operate with Vertex/Coordinate pairs, which enables a more
    efficient access to values and gradients by storing information
    directly in the vertex.

-   The concepts `InterpolationTraits` and `GradientFittingTraits`
    have been updated to reflect the real needs of the code (some
    types and operators were used in the code but did not appear in
    the concepts).

### CGAL and the Boost Graph Library (BGL)

-   Added a helper function, `CGAL::is_valid_polygon_mesh`, that
    checks the validity of a polygon mesh using BGL functions.

-   Improved the function `CGAL::Euler::collapse_edge` such that the
    target vertex of the collapsed edge is now always kept after the
    collapse.

-   The function `copy_face_graph()` now uses named parameters, some
    allowing it to use property maps instead of output iterators.

-   Addition of the following named parameters :
    -   vertex_to_vertex_output_iterator
    -   halfedge_to_halfedge_output_iterator
    -   face_to_face_output_iterator
    -   vertex_to_vertex_map
    -   halfedge_to_halfedge_map
    -   face_to_face_map

### CGAL and Solvers

-   **Breaking change**: `CGAL::Diagonalize_traits` is now deprecated
    and should not be used. The class `CGAL::Eigen_diagonalize_traits`
    (along with the _Eigen_ library) should be used instead.

### CGAL and Boost Property Maps

-   Added a read-write property map to convert on-the-fly geometric
    objects from Cartesian kernels.

### 2D Arrangements

-   Refracted and fixed the `graph_traits` for the dual of an arrangement of the
    following types:
    `Arrangement_on_surface_2`,
    `Arrangement_2`,
    `Arrangement_on_surface_with_history_2`, and
    `Arrangement_with_history_2`.

-   **Breaking change**: The old `<CGAL/graph_traits_Dual_Arrangement_2.h>`
    header file has been replaced by the four header files below; each defines
    the `graph_traits` for dual of the corresponding arrangement type.
    `<CGAL/graph_traits_dual_arrangement_on_surface_2.h>`,
    `<CGAL/graph_traits_dual_arrangement_2.h>`,
    `<CGAL/graph_traits_dual_arrangement_on_surface_with_history_2.h>`, and
    `<CGAL/graph_traits_dual_arrangement_with_history_2.h`.


Release 4.12
------------

Release date: April 2018


### Important Notice

-   The CMake scripts used by CGAL have been changed to use modern patterns
    introduced by CMake 2.8.12 and CMake 3.0: instead of setting CMake
    variables, the script now defines imported targets and uses link
    interfaces.

    That is mostly backward-compatible with existing usages of CGAL CMake
    scripts. The only non-compatible effect is that the `CMAKE_BUILD_TYPE`
    and compilation flags are no longer copied from the `CGAL_DIR` to the
    project using it. Note also that the `CMAKE_BUILD_TYPE` is no longer
    set to `Release` by default. For a developer using the Visual Studio
    IDE or the Xcode IDE, the change should be transparent. Developers using
    makefiles or the Ninja build-system should set the `CMAKE_BUILD_TYPE`
    to `Release` manually, to avoid using CGAL libraries without any
    compile-time optimization.

### Header-only Mode

-   Since CGAL-4.9, it has been possible to use CGAL by configuring it
    using CMake, but without compiling the CGAL libraries. With CGAL-4.12,
    it is now possible to use CGAL header-only, without even configuring
    it. CMake is then used only to configure programs using CGAL.

### Compiler Support

-   The Microsoft Visual C++ 2017 version 15.3 has introduced support for
    C++17, with the compilation flag `/std:c++17`. CGAL 4.12 has an initial
    support for that flag: the code will compile, but a lot of deprecation
    warnings will remain. Note that Boost version 1.67 is the first version
    of Boost supporting `/std:c++17`.

-   The compilation flag `/permissive-` of Visual C++ is now supported.

### 2D Movable Separability of Sets (new package)

-   A new package called "2D Movable Separability of Sets" has been
    introduced. It handles a class of problems that deal with moving
    sets of objects in the plane; the challenge is to avoid collisions
    between the objects while considering different kinds of motions and
    various definitions of separation.

    At this point this package consists of the implementations of
    various predicates and constructions related to castings of
    polygonal objects. In particular, it can be used to determine
    whether a feasible mold for a polygonal object does exist. If a mold
    exists, the package can also be used to compute all possible
    orientations of the feasible molds and the corresponding motions
    needed to remove the casted object from the mold.

### Classification (new package)

-   This package offers an algorithm that classifies a data set into a
    user-defined set of labels (such as ground, vegetation, buildings,
    etc.). A flexible API is provided so that users can classify any
    type of data, compute their own local features on the input data
    set, and define their own labels.

### Kinetic Data Structures (removed package)

-   This package has been removed from CGAL-4.12. Users of the package
    will have to keep using the source code available in CGAL-4.11 or
    earlier.


### 3D Convex Hull

-   **Breaking change**: The header `<CGAL/convex_hull_3.h>` no longer
    includes `<CGAL/Polyhedron_3.h>`, as the convex hull function works
    with any model of the concept `MutableFaceGraph`.

### 2D Arrangements

-   When removing an edge from an arrangement and the user has requested to
    remove the end-vertices in case they become redundant (either isolated or
    approach infinity), defer the removal of the such end-vertices to occur
    after the observer is notified that the edge has been removed. This is
    symmetric (opposite) to the order of notification when an edge is inserted.

    The user can restore old (non-symmetric) behavior by defining the macro:

    `CGAL_NON_SYMETRICAL_OBSERVER_EDGE_REMOVAL_BACKWARD_COMPATIBILITY`

### 2D Periodic Triangulations

-   **Breaking change**: The class
    `Periodic_2_triangulation_hierarchy_vertex_base_2` (and its
    corresponding header) have been removed. Users should directly use
    the class `Triangulation_hierarchy_vertex_base_2`, which is
    identical.
-   **Breaking change**: The functions `circumcenter()`,
    `side_of_oriented_circle()`, and `is_extensible_in_1_sheet_h[12]()`
    are related to Delaunay triangulations and have been moved from
    `Periodic_2_triangulation_2` to
    `Periodic_2_Delaunay_triangulation_2`.

### 2D Alpha Shapes

-   It is now possible to use `CGAL::Periodic_2_triangulation_2` as
    underlying triangulation for `Alpha_shape_2`.

### 3D Surface Mesh Generation

-   Add the function `facets_in_complex_2_to_triangle_mesh()` that
    exports `Surface_mesh_complex_2_in_triangulation_3` facets into
    a `MutableFaceGraph`.

### 3D Mesh Generation

-   Add the function `facets_in_complex_3_to_triangle_mesh()` that
    exports `Mesh_complex_3_in_triangulation_3` facets into a
    `MutableFaceGraph`.
-   **Breaking change:** The concept `MeshDomainWithFeatures_3` has been
    modified, to improve the performance and the reliability of the
    sampling of 1D curves of the domain.
-   Add the ability to ensure that the output mesh surface describes a
    manifold, when the input surface is a manifold. New named parameters
    `manifold()`, `manifold_with_boundary()`, and `non_manifold()` are
    added.

### Optimal Transportation Curve Reconstruction

-   New method `run_under_wasserstein_tolerance()` which allows the
    user to perform curve reconstruction by relying on a threshold on
    the Wasserstein distance. This is useful when the number of edges
    in the expected output reconstruction is not known.

### Polygon Mesh Processing

-   Added two functions for orienting connected components :
    -   `CGAL::Polygon_mesh_processing::orient()`
    -   `CGAL::Polygon_mesh_processing::orient_to_bound_a_volume()`

-   Added a new function for intersection tests between triangle meshes
    and/or polylines or range of polylines, and another one to report
    all the pairs of meshes intersecting from a range of meshes:
    -   `CGAL::Polygon_mesh_processing::do_intersect()`
    -   `CGAL::Polygon_mesh_processing::intersecting_meshes()`

-   Added new functions for feature detection and feature-guided
    segmentation:
    -   `CGAL::Polygon_mesh_processing::detect_sharp_edges()`
    -   `CGAL::Polygon_mesh_processing::detect_vertex_incident_patches()`
    -   `CGAL::Polygon_mesh_processing::sharp_edges_segmentation()`

### Point Set Shape Detection

-   **Breaking change**:
    `CGAL::Shape_detection_3::Efficient_RANSAC_traits` is now called
    `CGAL::Shape_detection_3::Shape_detection_traits`.
-   New algorithm: `CGAL::Region_growing`. This is a deterministic
    alternative to RANSAC for plane detection.
-   **Breaking change**: the API of `CGAL::regularize_planes()` is
    generalized to accept other types of input than the RANSAC output.
-   Add a callback mechanism for both `CGAL::Efficient_RANSAC` and
    `CGAL::Region_growing`.

### Point Set Processing

-   **Breaking change**: the API of `CGAL::structure_point_set()` is
    generalized to accept other types of input than the RANSAC output.
-   **Breaking change**: the API of all functions of Point Set
    Processing is modified to use ranges (instead of iterators) and
    Named Parameters (similarly to the API of Polygon Mesh
    Processing). The old API is kept as deprecated.

### CGAL and the Boost Graph Library (BGL)

-   Add helper function `CGAL::expand_face_selection_for_removal` that
    expands a face selection to avoid creating a non manifold mesh when
    removing the selected faces.

-   Added support for dynamic property maps.

-   Added an interface to the [METIS library], which allows to partition
    any mesh that is a model of `FaceListGraph`.  Wrappers to the
    METIS functions `METIS_PartMeshNodal` and `METIS_PartMeshDual` are
    offered.

    [METIS library]: http://glaros.dtc.umn.edu/gkhome/metis/metis/overview


Release 4.11
------------

Release date: September 2017

### 3D Periodic Regular Triangulations (new feature)

-   Added the class `Periodic_3_regular_triangulation_3`, which provides
    functionality for 3D periodic weighted Delaunay triangulations. The
    construction is fully dynamic: it provides both point insertion and
    vertex removal.

### dD Regular Triangulations (new feature)

-   Added the class `Regular_triangulation`, which provides
    functionality for dD weighted Delaunay triangulations. Note that the
    removal of points is not yet supported.

### 2D and 3D Linear Geometry Kernel (breaking change)

-   **Breaking change**: The dangerous implicit conversions between
    weighted points and points in the concept `Kernel` have been
    disabled. Constructors offering to build a weighted point from a
    point (and reversely) are still requested by the concept `Kernel`
    but must now be marked with the `explicit` specifier.
-   **Breaking change**: The removal of implicit conversions between
    points and weighted points in the concept `Kernel` has incidentally
    created various minor breaking changes in the following packages: 2D
    Alpha Shapes, 2D and 3D Triangulations, and 3D Mesh Generation. See
    the full changelog for details.

### Surface Mesh

-   **Breaking change**:
    `operator >>(std::istream&,       Surface_mesh&)` no longer clears
    the surface mesh.

### Triangulated Surface Mesh Parameterization (breaking change)

-   **Breaking change**: The package has been rewritten and can operate
    on any model of the `MutableFaceGraph` concept. All previous
    parameterization methods are still offered, although with a
    different, simpler API. The documentation has been updated and
    offers a gentle introduction to the new API. Users who wish to use
    the former API must use a version prior to 4.11.
-   **Breaking change**: The adapter to add virtual seams is now the
    class `CGAL::Seam_mesh` in the package *CGAL and the BGL*.
-   **Breaking change**: The package has been restructured and most
    headers have been moved. In a general manner, users should replace
    `<CGAL/XXX.h>` with `<CGAL/Surface_mesh_parameterization/XXX.h>`
-   Add the *As Rigid As Possible Parameterization* method. This
    parameterization allows the user to prioritize angle preservation,
    shape preservation, or a balance of both.
-   Add the *Orbifold Tutte Embedding* method. This parameterization
    method allows to parameterize meshes that are topological spheres.

### 3D Surface Subdivision Methods (breaking changes)

-   The subdivision algorithms now work on any model of a
    `MutableFaceGraph`. A new API to the subdivision methods is offered,
    which uses optional named parameters to pass the number of
    iterations and a vertex property map.
-   **Breaking change**: Removed the headers
    `<CGAL/Subdivision_method_3.h>` and `<CGAL/Subdivision_mask_3.h>`.
    The headers `<CGAL/Subdivision_method_3/subdivision_methods_3.h>`
    and `<CGAL/Subdivision_method_3/subdivision_masks_3.h>` should
    respectively be used instead.
-   Sqrt3 subdivision can now handle input surfaces with a border.

### Scale-Space Surface Reconstruction (breaking change)

-   **Breaking change**: the API was rewritten to separate the smoothing
    and meshing algorithm and making it possible for the user to use
    different ones. The default algorithms used are the same as before
    this API change, but methods are moved to the classes
    `Weighted_PCA_smoother` and `Alpha_shape_mesher`.
-   Alternative smoothing and meshing methods are provided:
    `Jet_smoother` and `Advancing_front_mesher`.

### 2D Alpha Shapes

-   **Breaking change**: Mirrored the concepts of the 2D alpha shape
    package with those of the 3D Alpha Shapes package. Consequently, a
    new concept, `WeightedAlphaShapeTraits_2`, is introduced to provide
    requirements on the traits class for 2D weighted alpha shapes. All
    models of the concept `Kernel` are models of this new concept.
-   The concept `AlphaShapeTraits_2` now provides requirements on the
    traits class for 2D basic alpha shapes, and refines
    `DelaunayTriangulationTraits_2`.

### Interpolation

-   **Breaking change**: The concept `GradientFittingTraits` now
    additionally requests a weighted point type `Weighted_point_d` and a
    functor `Construct_point_d`. The model
    `CGAL::Interpolation_gradient_fitting_traits_2` has been
    appropriately modified to still be a model of the concept
    `GradientFittingTraits`.

### 2D and 3D Triangulations

-   **Breaking change**: Added a new functor requirement,
    `Construct_point_2`, to the concepts `TriangulationTraits_2` and
    `RegularTriangulationTraits_2` and a new functor requirement,
    `Construct_point_3`, to the concepts `TriangulationTraits_3` and
    `RegularTriangulationTraits_3`. All models of the concept `Kernel`
    already provide these functors.
-   **Breaking change**: Introduced the concepts
    `RegularTriangulationVertexBase_2` and
    `RegularTriangulationVertexBase_3`. These concepts describe the
    requirements on classes meant to represent a vertex of a regular
    triangulation. Concepts that previously refined
    `TriangulationVertexBase_2` or `TriangulationVertexBase_3` but
    described in fact a vertex class used in a regular triangulation,
    such as the concept `MeshVertexBase_3` in the 3D mesh generation
    package, now refine the corresponding new regular vertex concept.
-   **Breaking change**: Uniformized the point type across all vertex
    and cell concepts. The triangulation point type name is now always
    `Point`. Note that this does not change the requirements but only
    the name: `Point` is still expected to be equal to
    `Traits::Point_[23]` for basic and Delaunay triangulations or to
    `Traits::Weighted_point_[23]` for regular triangulations.
    Consequently:
    -   The concept `RegularTriangulationVertexBase_2` now requests a
        `Point` type (equal to `Traits::Weighted_point_2`)
    -   The concept `RegularTriangulationCellBase_3` now requests a
        `Point` type instead of a `Weighted_point` type (but still equal
        to `Traits::Weighted_point_3`)
    -   The concept `DelaunayTriangulationCellBase_3` now requests a
        `Point` type instead of a `Point_3` type (but still equal to
        `Traits::Point_3`).
-   Introduced a new concept,
    `RegularTriangulationCellBaseWithWeightedCircumcenter_3`, which
    describes the requirements on a cell of a regular triangulation that
    caches its circumcenter. The existing class
    `Regular_triangulation_cell_base_with_weighted_circumcenter_3` is
    the default model of this concept.
-   Added a new 3D traits class,
    `Robust_weighted_circumcenter_filtered_traits_3` which provides
    robust versions of the kernel functors
    `Construct_weighted_circumcenter_3`, `Compute_squared_radius_3`, and
    `Compute_squared_radius_smallest_orthogonal_sphere_3`. This class
    can be used as traits class in the `Mesh_3` package to
    efficiently yet robustly generate 3D meshes.
-   Add a new type of polyhedral domain with features,
    `Polyhedral_complex_mesh_domain_3`. The domain is defined by a
    collection of triangulated surfaces, forming a complex.

### 3D Periodic Triangulations

-   Added new locate and geometric access functions for 3D periodic
    triangulations.
-   The class `Periodic_3_Delaunay_triangulation_traits_3` now inherits
    `Periodic_3_triangulation_traits_3`.
-   **Breaking change**: Some geometric access functions in
    `Periodic_3_triangulation_3` were renamed. The introduction of
    `Periodic_3_regular_triangulation_3` required to distinguish between
    functions such as `segment()` returning a segment of weightless
    points, or a segment of weighted points. As a general rule, previous
    geometrical access functions will return objects with point type
    that of the triangulation (thus, weighted objects when using
    weighted triangulations) and functions containing `construct` in the
    name will always return weightless geometrical objects.
-   **Breaking change**: The concept `Periodic_3TriangulationTraits_3`
    now requests a domain getter: `get_domain()`.
-   Introduced a new concept,
    `RegularTriangulationCellBaseWithWeightedCircumcenter_3`, which
    describes the requirements on a cell of a regular triangulation that
    caches its circumcenter. The existing class
    `Regular_triangulation_cell_base_with_weighted_circumcenter_3` is
    the default model of this concept.

### 3D Mesh Generation

-   **Breaking change**: The type of the surface center in the concept
    `MeshCellBase_3` has been changed from `Triangulation::Point` to
    `TriangulationTraits::Point_3` to reflect that it is a weightless
    point.
-   **Breaking change**: The function `invalidate_circumcenter()` of the
    concept `MeshCellBase_3` is renamed to
    `invalidate_weighted_circumcenter_cache()` and moved to the new
    concept `RegularTriangulationCellBaseWithWeightedCircumcenter_3`,
    which the concept `MeshCellBase_3` now refines.

### Poisson Surface Reconstruction

-   A new global function
    `CGAL::poisson_surface_reconstruction_delaunay()` is provided in
    addition to the current class-based API in order to make it easier
    to use.

### Point Set Processing

-   New functions to read from and write to LAS/LAZ files (LIDAR
    format), with or without taking additional properties into account.
-   **Breaking change:** The API of the PLY function to read points with
    properties is modified for unification with LAS (see
    `CGAL::read_ply_points_with_properties()`). A new function to write
    PLY with properties is provided
    (`CGAL::write_ply_points_with_properties()`).

### Spatial Searching

-   Add function `Kd_tree::remove(Point)`.

### 3D Fast Intersection and Distance Computation

-   Add a template parameter to `AABB_traits` for a property map that
    associates a bounding box to a primitive

### CGAL and the Boost Graph Library

-   Add a partial specialization for the class
    `CGAL::Linear_cell_complex_for_combinatorial_map` so that it is a
    model of the graph concepts `BidirectionalGraph` and
    `EdgeAndVertexListGraph` and of the concept `MutableFaceGraph`. This
    class can thus now be used in all BGL functions and algorithms.
-   Helper functions to create an icosahedron, a regular prism and a
    pyramid have been added.
-   Add class `CGAL::Face_filtered_graph` that wraps an existing graph
    and hide all simplices that are not in the selected connected
    components.
-   Added the class `CGAL::Seam_mesh`. The `Seam_mesh` is a graph
    adaptor which allows to create virtual borders when marking edges as
    seam edges.
-   Add the functions `read_off()` and `write_off()`.

Release 4.10
------------

Release date: May 2017

### Installation

-   The minimum required version of CMake is now 3.1. All CMake versions
    up to 3.7 are supported.
-   The compilation of some demo may require a C++11 compiler. The CGAL
    library itself still support C++03 compilers.
-   The shell script `cgal_create_cmake_script` now enables C++14 by
    default.
-   A new mechanism to check which packages of CGAL are used have been
    added. It is particularly interesting for commercial users to ensure
    they have a valid commercial license for the packages they used.
    This can also be used to make sure only LGPL header files are used.
-   Because of a bug in the g++ compiler about the C++11 keyword
    `thread_local`, the CGAL\_Core library now always requires
    `Boost.Thread` if the g++ compiler is used.

### Generalized Maps (new package)

-   This package implements Generalized Maps in d dimensions. A
    generalized map is a data structure enabling to represent an
    orientable or non orientable subdivided object by describing all the
    cells of the subdivision (for example in 3D vertices, edges, faces,
    volumes) and all the incidence and adjacency relationships between
    these cells. This data structure is the generalization of the
    combinatorial maps in order to be able to represent non orientable
    objects.

### 3D Point Set (new package)

-   This package provides a flexible data structure `CGAL::Point_set_3`
    that allows the user to easily handle point sets with an arbitrary
    number of attributes (such as normal vectors, colors, labeling,
    etc.).

### Combinatorial Maps and Linear cell complex

-   **Breaking change**: the requirements of the item class used to
    customize a combinatorial map and a linear cell complex changed.
    Instead of defining the type of darts used, you have to define the
    information you want to add in each dart. You can define the
    `CGAL_CMAP_DART_DEPRECATED` macro to keep the old behavior.

### Triangulated Surface Mesh Shortest Paths

-   **Breaking change**: Rename all functions, types, and enums using
    *barycentric coordinate* to *barycentric coordinates*.

### CGAL and the Boost Graph Library (BGL)

-   **Breaking change**: Addition of a free function `reserve()` in the
    concept `MutableFaceGraph`. Models provided by CGAL have been
    updated.

### 2D and 3D Linear Geometry Kernel

-   **Breaking change**: The function `compare_slopes()` was renamed
    `compare_slope`.
-   Added a 2D and 3D weighted point class and predicates and
    constructions.
-   Add functions `l_infinity_distance()` for 2D and 3D.
-   Add a new functor in CGAL Kernel concept to compare the slope of two
    3D segments. All models of the Kernel concept now provide the
    functor `Compare_slope_3`, and the free function `compare_slope()`
    is available.
-   Add an operator in CGAL Kernel concept `Angle_3` to qualify the
    angle between the normal of the triangle given by three points, and
    a vector.

### 3D Convex Hull

-   The convex hull function can also produce a `Surface_mesh`, and
    generally speaking any model of the concept `MutableFaceGraph`
-   The function `convex_hull_3_to_polyhedron_3()` is deprecated and
    `convex_hull_3_to_face_graph.h` should be used instead.
-   The class `Convex_hull_traits_3` now documents a nested type
    `Polygon_mesh` instead of `Polyhedron_3`. The other nested type is
    kept for backward compatibility.
-   Remove the function `CGAL::convex_hull_incremental_3()` deprecated
    since CGAL 4.6.

### 3D Boolean Operations on Nef Polyhedra

-   Add a new constructor from a face graph model

### Linear cell complex

-   Deprecate class `Linear_cell_complex` which is now renamed
    `Linear_cell_complex_for_combinatorial_map_dart`.

### 2D Triangulation data structure

-   Add function `insert_in_hole`.

### 2D Triangulations

-   **Breaking change**: Removed the arbitrary dimensional weighted
    point class. Users must use a version prior to 4.9 if they need this
    class.
-   **Breaking change**:The number type of weighted points in regular
    triangulations is no longer a template parameter but the field type
    of the geometric traits class. Users who need this feature must use
    a version prior to 4.9
-   The class `Regular_triangulation_filtered_traits_2` deprecated since
    CGAL 3.6 has been removed.
-   Deprecate the class `Regular_triangulation_euclidean_traits_2`, as
    the weighted point and the function objects for weighted points are
    part of the concept `Kernel`/
-   The class `Regular_triangulation_2` can take a kernel as template
    argument, that is one no longer has to use
    `Regular_triangulation_euclidea_traits_2`, although this still
    works.

### 3D Triangulations

-   **Breaking change**: The number type of weighted points in regular
    triangulations is no longer a template parameter but the field type
    of the geometric traits class. Users who need this feature must use
    a version prior to 4.9.
-   The class `Regular_triangulation_filtered_traits_3` deprecated since
    CGAL 3.6 has been removed.
-   Deprecate the class `Regular_triangulation_euclidean_traits_3`, as
    the weighted point and the function objects for weighted points are
    part of the concept `Kernel`/
-   The class `Regular_triangulation_3` can take a kernel as template
    argument, that is one no longer has to use
    `Regular_triangulation_euclidean_traits_3`, although this still
    works.
-   Add function `link_to_face_graph()` to copy the set of faces
    incident to a vertex into a model of `FaceGraph`.

### 3D Mesh Generation

-   The constructor
    `CGAL::Polyhedral_mesh_domain_with_features_3(std::string)` is
    deprecated.

### Polygon Mesh Processing

-   Add fast and robust corefinement and Boolean operation functions for
    triangulated surface meshes:
    -   `CGAL::Polygon_mesh_processing::corefine_and_compute_union()`
    -   `CGAL::Polygon_mesh_processing::corefine_and_compute_difference()`
    -   `CGAL::Polygon_mesh_processing::corefine_and_compute_intersection()`
    -   `CGAL::Polygon_mesh_processing::corefine()`
    -   `CGAL::Polygon_mesh_processing::does_bound_a_volume()`
    -   `CGAL::Polygon_mesh_processing::surface_intersection()`
-   Add functions to compute approximated Hausdorff distances between
    two meshes, a mesh and a point set, or a point set and a mesh:
    `sample_triangle_mesh()`, `approximated_Hausdorff_distance()`,
    `approximated_symmetric_Hausdorff_distance()`,
    `max_distance_to_triangle_mesh()`, `max_distance_to_point_set()`.
-   The function `CGAL::Polygon_mesh_processing::bbox_3()` has been
    renamed `CGAL::Polygon_mesh_processing::bbox()`.

### Point Set Processing

-   Function `CGAL::remove_outliers()` has an additional parameter based
    on a distance threshold to make it easier and more intuitive to use.
-   New functions for automatic scale estimations: either a global scale
    or multiple local scales can be estimated for both 2D and 3D point
    sets based on the assumption that they sample a curve in 2D or a
    surface in 3D.

### CGAL and the Boost Graph Library (BGL)

-   Add function `CGAL::convert_nef_polyhedron_to_polygon_mesh()` to
    convert a `Nef_polyhedron_3` to any model of the `MutableFaceGraph`
    concept.
-   Add class `CGAL::Graph_with_descriptor_with_graph` that wraps an
    existing graph and provides a reference to the said graph to all of
    its descriptors.

### Cone Based Spanners

-   Add a parameter to compute half Tao graph and half Theta graph.
-   Add an ipelet for this package.

### Geometric Object Generators

-   Add point random generators
    -   in a 3D triangle mesh model of the concept `FaceListGraph`
        (`CGAL::Random_points_in_triangle_mesh_3`),
    -   on the boundary of a tetrahedral mesh
        (`CGAL::Random_points_in_tetrahedral_mesh_boundary_3`),
    -   in a tetrahedral mesh
        (`CGAL::Random_points_in_tetrahedral_mesh_3`),
    -   in a 2D triangle mesh
        (`CGAL::Random_points_in_triangle_mesh_2`),
    -   in a range of 2D or 3D triangles
        (`CGAL::Random_points_in_triangles_3` and
        `CGAL::Random_points_in_triangles_2`).
    -   on a 3D segment (`CGAL::Random_points_on_segment_3`).

Release 4.9
-----------

Release date: Sept 2016

### Header-only mode

-   CGAL can now be used in headers only mode, i.e. without compiling
    the CGAL libraries and linking with these libraries when compiling
    examples, tests and demos. Note that running CMake on CGAL is still
    required in order to generate some configuration files.

### Cone Based Spanners (new package)

-   This package provides algorithms for constructing two kinds of
    cone-based spanners: Yao graph and Theta graph, given a set of
    vertices on the plane and the directions of cone boundaries.

### 2D Minkowski Sums

-   Introduce a convex decomposition strategy, namely
    `Polygon_nop_decomposition_2`, that merely passed the input polygon
    to the list of output polygons.
-   Introduce overloads of the function `minkowski_sum_2()`, which
    accepts 2 decomposition strategies.
-   Introduce an overloaded function called
    `minkowski_sum_by_decomposition_2(P, Q, decom_no_holes,     decomp_with_holes)`,
    which computes the 2D Minkowski sum using optimal choices of
    decomposition strategies.

### Combinatorial Maps

-   Deprecate global functions (`make_combinatorial_hexahedron()`,
    `make_combinatorial_polygon()`, `make_combinatorial_tetrahedron()`,
    `make_edge()`, `insert_cell_0_in_cell_1()`,
    `insert_cell_0_in_cell_2()`, `insert_cell_1_in_cell_2()`,
    `insert_cell_2_in_cell_3()`, `insert_dangling_cell_1_in_cell_2()`,
    `is_insertable_cell_1_in_cell_2()`,
    `is_insertable_cell_2_in_cell_3()`, `is_removable()`,
    `remove_cell()`) which are now member functions in the
    `CombinatorialMap` concept.
-   It is not longer possible to use the old API switched on by defining
    the macro `CGAL_CMAP_DEPRECATED`. This API was deprecated since CGAL
    4.4.

### Point Set Processing

-   New function `CGAL::read_ply_custom_points()` that allows the user
    to read any additional point attribute from a PLY input point set.
-   `CGAL::structure_point_set()`: new algorithm that takes advantage of
    detected planes to produce a structured point set (with flat
    regions, sharp edges and vertices).

### Point Set Shape Detection

-   New post-processing algorithm: `CGAL::regularize_planes()`. This
    allows the user to favor parallelism, orthogonality, coplanarity
    and/or axial symmetry between detected planes.

### Polygon Mesh Processing

-   Add the function
    `CGAL::Polygon_mesh_processing::is_polygon_soup_a_polygon_mesh()` to
    check whether a polygon soup is a polygon mesh.
-   Add some new features to `CGAL::isotropic_remeshing()`:
    -   It is now possible to select fixed vertices that survive the
        remeshing process, and to keep face attributes such as colors
        valid after remeshing.
    -   The user can choose the number of relaxation steps happening at
        each loop, and to run 1-dimensional relaxation along constrained
        polylines.
-   The functions `CGAL::Polygon_mesh_processing::triangulate_face()`
    and `CGAL::Polygon_mesh_processing::triangulate_faces()` now
    indicate whether some faces have not been triangulated.

### Surface Mesh Deformation

-   Add a new tag `SRE_ARAP` to use the Smoothed Rotation Enhanced
    As-Rigid-As-Possible deformation algorithm.

### 3D Fast Intersection and Distance Computation

-   Add the functions `AABB_tree::first_intersection()` and
    `AABB_tree::first_intersected_primitive()` that compute the
    intersection which is closest to the source of a ray

### CGAL and the Boost Graph Library (BGL)

-   Add a helper function `CGAL::copy_face_graph()` to copy a source
    FaceListGraph into another FaceListGraph of different type.
-   Add a class `CGAL::Dual` that creates the dual view of a `FaceGraph`
    and a creation function `CGAL::dual(primal)`.

#### CGAL and Boost Property Maps

-   It is not longer possible to use the old API of the property maps
    provided by CGAL, switched on by defining the macro
    `CGAL_USE_PROPERTY_MAPS_API_V1`. This API was deprecated since CGAL
    4.3.

Release 4.8
-----------

Release date: April 2016

### General

-   The support for Qt3 is dropped and all demos using it got removed.

### Installation

-   Starting with Visual C++ 2015 we no longer require `Boost.Thread` as
    we use the C++11 keyword `thread_local` and the C+11 class
    `std::mutex` .
-   The same holds for g++ 4.8 or later when the C++11 standard is used.

### Optimal Transportation Curve Reconstruction (new package)

-   This package implements a method to reconstruct and simplify 2D
    point sets. The input is a set of 2D points with mass attributes,
    possibly hampered by noise and outliers. The output is a set of line
    segments and isolated points which approximate the input points.

### 2D Regularized Boolean Set-Operations

-   Improve the performance of operations in some settings.
    **Breaking change**: This improvement requires changes of the face
    and halfedge type of the underlying arrangement Dcel. See the
    concepts `GeneralPolygonSetDcelHalfedge` and
    `GeneralPolygonSetDcelFace` for more details. If you use a different
    simplex types, inheriting your simplices from `CGAL::Gps_face_base`
    and `CGAL::Gps_halfedge_base` is sufficient to accommodate for the
    update.

### 3D Boolean Operations on Nef Polyhedra

-   Add 3 new constructors: from a point range, from a point, and from a
    segment.

### Combinatorial Maps

-   **Breaking change**: Change the type of Boolean marks, old type is
    int, new type is `size_type`. If no more mark is available,
    `get_new_mark` throws an exception, instead of returning `-1`.

### 2D Arrangements

-   Speed up the edge removal in case the incident faces contains many
    holes.
-   Set the format of polylines and polycurves. The format of a general
    polyline or polycurve consists of the sequence of subcurves that
    comprise the original curve. The format of a polyline of linear
    segments consists of the sequence of points that define the original
    curve. (The latter restores the format used before polycurves were
    introduced in 4.7.) Fix the extraction from istream and insertion
    into ostream operators of polylines and polycurves accordingly.
-   Fix the traits class that handles Bezier curves. In particular, fix
    the case where the curve is closed (i.e, the first and last control
    points coincide).

### 3D Mesh Generation

-   Add support of 3D gray level images as input for the tetrahedral
    mesh generation.
-   **Breaking change:** All models of the concept `MeshDomain_3` must
    now provide a member function `bbox()`.

### Advancing Front Surface Reconstruction

-   Optional template functor `Filter` is replaced by another optional
    template functor `Priority`. This allows to change the way facets
    are prioritized by the algorithm instead of having a simple option
    to reject some facets.
    **Breaking change**: Programs using the old `Filter` API will not
    compile anymore as it must be replaced with the `Priority` API as
    described in the manual. Codes using the default behavior are not
    impacted.

### Polygon Mesh Processing

-   Add a new triangle-based isotropic remeshing algorithm for
    triangulated surface meshes,
    `CGAL::Polygon_mesh_processing::isotropic_remeshing()` and a helper
    function for isotropic remeshing :
    `CGAL::Polygon_mesh_processing::split_long_edges()`
-   Add the function `CGAL::Polygon_mesh_processing::border_halfedges()`
    to collect the border of a given face range
-   Add the function
    `CGAL::Polygon_mesh_processing::remove_isolated_vertices()` to be
    used on any polygon mesh
-   Add the function `CGAL::Polygon_mesh_processing::triangulate_face()`
    to triangulate a single face of a polygon mesh
-   Add an overload for
    `CGAL::Polygon_mesh_processing::triangulate_faces()` to triangulate
    a range of faces of a polygon mesh
-   Add function `keep_large_connected_components()`
-   Add measuring functions for polygon meshes, to compute length, area,
    and volume of simplices or group of simplices of a polygon mesh.
-   Add function `bbox_3()` to compute the bounding box of a polygon
    mesh.

### Point Set Processing

-   **Breaking change:** new template parameter `Concurrency_tag` for
    the functions `compute_average_spacing()`,
    `edge_aware_upsample_point_set()`, `jet_estimate_normals()`,
    `jet_smooth_point_set()`, and `pca_estimate_normals()`. To update
    your code simply add as first template parameter
    `CGAL::Sequential_tag` or `CGAL::Parallel_tag` when calling one of
    these functions.
-   `CGAL::Parallel_tag` can no longer be used in Point Set Processing
    algorithms if TBB is not available.
-   Add a new simplification algorithm based on hierarchical clustering:
    `CGAL::hierarchy_simplify_point_set()`. It allows either to
    uniformly simplify the point set or to automatically adapt the local
    density of points to the local variation of the input computed by
    principal component analysis.
-   New IO functions for PLY format (Polygon File Format):
    `CGAL::read_ply_points()`, `CGAL::read_ply_points_and_normals()`,
    `CGAL::write_ply_points()` and
    `CGAL::write_ply_points_and_normals()`.

### Surface Mesh Parameterization

-   `LSCM_parameterizer_3` now uses by default Eigen instead of OpenNL
    as a model of `SparseLinearAlgebraTraits_d`.

### Spatial Searching

-   Add function to find any point in a range query, that is neither all
    points, nor the closest one.

### Principal Component Analysis

-   Add a template parameter `DiagonalizeTraits` for functions
    `CGAL::linear_least_squares_fitting_2()` and
    `CGAL::linear_least_squares_fitting_3()`. This allows to either
    choose the legacy internal diagonalization code from CGAL or the
    Eigen implementation (or any class that is a model of
    `DiagonalizeTraits`). Variants of the function that automatically
    deduce the kernel also automatically select the diagonalizer, so the
    API is mostly preserved.

### CGAL and Solvers

-   This package now includes all CGAL concepts for solvers with models
    using the third party Eigen library.

### CGAL and the Boost Graph Library (BGL)

-   Add function `CGAL::split_graph_into_polylines()` that allows to
    extract from a soup of segments given as a graph, polylines with
    nodes of degree at most 2. In addition a functor can be passed to
    the function to specify additional polyline endpoints.
-   New functions to manipulate selection of faces, edges and vertices
    in a halfedge graph are added: `CGAL::expand_face_selection()`,
    `CGAL::reduce_face_selection()`, `CGAL::expand_edge_selection()`,
    `CGAL::reduce_edge_selection()` `CGAL::expand_vertex_selection()`,
    `CGAL::reduce_vertex_selection()` and
    `CGAL::select_incident_faces()`.
-   Add a helper function `CGAL::clear` which clears a MutableFaceGraph
    efficiently and generically.

Release 4.7
-----------

Release date: October 2015

### Installation

-   The minimum required version of CMake is now 2.8.11. CMake versions
    3.1, 3.2, and 3.3 are supported.
-   All Qt4 demos have been updated and now require Qt5 to be compiled.
    Qt5 version 5.3 or higher is required. The support for Qt4 is
    dropped. To compile libCGAL\_Qt5 and demos, you must set the cmake
    or environment variable `Qt5_DIR` to point to the path to the
    directory containing the file `Qt5Config.cmake` created by your Qt5
    installation. If you are using the open source edition it should be
    `/path-to/qt-everywhere-opensource-src-<version>/qtbase/lib/cmake/Qt5`.
-   The code of the 3D demos now uses modern OpenGL, with shaders,
    instead of the fixed pipeline API of OpenGL-1.
-   The Microsoft Windows Visual C++ compiler 2015 (VC14) is now
    supported. However, since this compiler is not officially supported
    by Intel TBB 4.4 and Qt 5.5 (the latest versions available at the
    time of this release), the parallelism features of CGAL and Qt5
    demos will not work.

### L Infinity Segment Delaunay Graphs (new package)

-   The package provides the geometric traits for constructing the
    segment Delaunay graph in the max-norm (L Infinity). The traits also
    contain methods to draw the edges of the dual of the segment
    Delaunay graph in the max-norm i.e., the segment Voronoi diagram in
    the max-norm. The algorithm and traits rely on the segment Delaunay
    graph algorithm and traits under the Euclidean distance. The segment
    Voronoi diagram in the max-norm has applications in VLSI CAD.

### Advancing Front Surface Reconstruction (new package)

-   This package provides a greedy algorithm for surface reconstruction
    from an unorganized point set. Starting from a seed facet, a
    piecewise linear surface is grown by adding Delaunay triangles one
    by one. The most plausible triangles are added first, in a way that
    avoids the appearance of topological singularities.

### Triangulated Surface Mesh Shortest Paths (new package)

-   The package provides methods for computing shortest path on
    triangulated surface meshes. Given a set of source points on the
    surface, this package provides a data structure that can efficiently
    provides the shortest path from any point on the surface to the
    sources points. There is no restriction on the genus or the number
    of connected components of the mesh.

### Triangulated Surface Mesh Skeletonization (new package)

-   This package provides a (1D) curve skeleton extraction algorithm for
    a triangulated polygonal mesh without borders based on the mean
    curvature flow. The particularity of this skeleton is that it
    captures the topology of the input. For each skeleton vertex one can
    obtain its location and its corresponding vertices from the input
    mesh. The code is generic and works with any model of the
    \`FaceListGraph\` concept.

### 3D Point-Set Shape Detection (new package)

-   This package implements the efficient RANSAC method for shape
    detection, contributed by Schnabel et al. From an unstructured point
    set with unoriented normals, the algorithm detects a set of shapes.
    Five types of primitive shapes are provided by this package: plane,
    sphere, cylinder, cone and torus. Detecting other types of shapes is
    possible by implementing a class derived from a base shape.

### 2D Visibility (new package)

-   This package provides several variants to compute the visibility
    area of a point within polygonal regions in two dimensions.

### Polygon Mesh Processing (new package)

-   This package implements a collection of methods and classes for
    polygon mesh processing, ranging from basic operations on simplices,
    to complex geometry processing algorithms. The implementation of
    this package mainly follows algorithms and references given in
    Botsch et al.'s book on polygon mesh processing.

### General

-   Support for unordered sets and maps of the stdlib and of boost for
    handle and index classes.

### Approximation of Ridges and Umbilics on Triangulated Surface Meshes

-   This package now supports any model of the concept `FaceGraph`.
-   **Breaking change:** The package no longer supports models of
    `TriangulatedSurfaceMesh` which are not at the same time models of
    the concept `FaceGraph`.

### dD Geometry Kernel

-   Epick\_d gains 3 new functors: `Construct_circumcenter_d`,
    `Compute_squared_radius_d`, `Side_of_bounded_sphere_d`. Those are
    essential for the computation of alpha-shapes.

### 2D Arrangements

-   Introduced a new traits class, called
    `Arr_polycurve_traits_2<SubcurveTraits>`, which handles general
    piece-wise (polycurve) curves. The pieces do not necessarily have to
    be linear.
-   Introduced two new concepts called `ArrangementApproximateTraits_2`
    and `ArrangementConstructXMonotoneCurveTraits_2`.
    -   The existing `ArrangementLandmarkTraits_2` concept, which has
        two requirements, now refines the two respective concepts above.
    -   The template parameter of the existing
        `Arr_polyline_traits_2<SegmentTraits>` template must be
        substituted with a traits class that is a model of the
        `ArrangementConstructXMonotoneTraits_2` concept among the other
        when `Arr_polyline_traits_2` is instantiated.

### 2D Minkowski Sums

-   Added support for polygons with holes and optimized the construction
    of Minkowski sums.
    -   Introduced an implementation of the "reduced convolution"
        method, a variant of the method described in "2D Minkowski Sum
        of Polygons Using Reduced Convolution" by Behar and Lien. The
        new method supports polygons with holes and in many cases out
        performs the implementation of the existing (full) convolution
        method.
    -   Introduced two new classes that decompose polygons into convex
        pieces (models of the `PolygonConvexDecomposition_2` concept)
        based on vertical decomposition and constrained Delaunay
        triangulation, respectively. These new models also support the
        convex decomposition of polygons with holes.

### 3D Periodic Triangulations

-   Rename `Periodic_3_triangulation_traits_3` to
    `Periodic_3_Delaunay_triangulation_traits_3`.
-   Rename the concept `Periodic_3TriangulationTraits_3` to
    `Periodic_3DelaunayTriangulationTraits_3`.
-   Create `Periodic_3_triangulation_traits_3` and the concept
    `Periodic_3TriangulationTraits_3`.

### 2D Conforming Triangulations and Meshes

-   Add an optimization method `CGAL::lloyd_optimize_mesh_2()` that
    implements the Lloyd (or Centroidal Voronoi Tessellation)
    optimization algorithm in a Constrained Delaunay Triangulation. For
    optimization, the triangulation data structure on which the mesher
    relies needs its `VertexBase` template parameter to be a model of
    the new concept `DelaunayMeshVertexBase_2`.

### Point Set Processing and Surface Reconstruction from Point Sets

-   Add the function `CGAL::compute_vcm()` for computing the Voronoi
    Covariance Measure (VCM) of a point set. The output of this function
    can be used with the function `CGAL::vcm_is_on_feature_edge()` to
    determine whether a point is on or close to a feature edge. The
    former function is also internally used by
    `CGAL::vcm_estimate_normals()` to estimate the normals of a point
    set and it is particularly suited to point sets with noise.

### Spatial Sorting

-   Add the possibility to sort points on a sphere along a space-filling
    curve using the functions `CGAL::hilbert_sort_on_sphere` and
    `CGAL::spatial_sort_on_sphere`.

### Geometric Object Generators

-   Add new random generator of points in a 2D and 3D triangle and in a
    tetrahedron (`CGAL::Random_points_in_triangle_2`,
    `CGAL::Random_points_in_triangle_3`,
    `CGAL::Random_points_in_tetrahedron_3`).

Release 4.6.2
-------------

Release date: August 2015

This release only fixes bugs. See the list of fixed bugs on Github:

<https://github.com/CGAL/cgal/issues?q=milestone%3A4.6.2>

Release 4.6.1
-------------

Release date: June 2015

This release only fixes bugs. See the list of fixed bugs on Github:

<https://github.com/CGAL/cgal/issues?q=milestone%3A4.6.1+-label%3Ainvalid>

Release 4.6
-----------

Release date: April 2015

### Installation

-   The required version of Boost is now 1.48 or higher.

### 2D Polyline Simplification (new package)

-   This package enables to simplify polylines with the guarantee that
    the topology of the polylines does not change. This can be done for
    a single polyline as well as for a set of polyline constraints in a
    constrained triangulation. The simplification can be controlled with
    cost and stop functions.

### 2D Generalized Barycentric Coordinates (new package)

-   This package offers an efficient and robust implementation of
    two-dimensional closed-form generalized barycentric coordinates
    defined for simple two-dimensional polygons.

### Scale-Space Surface Reconstruction (new package)

-   This new package provides a class gathering a dedicated smoothing
    algorithm and some convenience functions to help the creation of a
    surface out of a point set using the 3D Alpha Shapes package. The
    particularity of this reconstruction pipeline is that the input
    point are in the output and no new points are created. Note that in
    the current version, the output is a triangle soup that is not
    necessarily a valid (manifold) polyhedral surface.

### Surface Mesh (new package)

-   The surface mesh class provided by this package is an implementation
    of the halfedge data structure allowing to represent polyhedral
    surfaces. It is an alternative to the packages `CGAL::Polyhedron_3`
    and `CGAL::HalfedgeDS`.

### dD Triangulation (new package)

-   This new package provides classes for manipulating triangulations in
    Euclidean spaces whose dimension can be specified at compile-time or
    at run-time. It also provides a class that represents Delaunay
    triangulations.

### dD Convex Hulls and Delaunay Triangulations

-   This package is deprecated and the new package Triangulation should
    be used instead.

### dD Geometry Kernel

-   It has been reported that the recently introduced `Epick_d` kernel
    may not work with Intel C++ Compiler prior to version 15.
    Documentation has been updated.

### 3D Convex Hulls

-   Add functions `halfspace_intersection_3` and
    `halfspace_intersection_with_constructions_3` to compute the
    intersection of halfspaces defining a closed polyhedron.
-   Fix a bug introduced in CGAL 4.5 that can appear while computing the
    convex hull of coplanar points.
-   Fix a robustness issue in `Convex_hull_traits_3`. This traits is
    used by default with the kernel
    `Exact_predicates_inexact_constructions_kernel`.
-   The function `CGAL::convex_hull_incremental_3` is deprecated and the
    function `convex_hull_3` should be used instead.

### Combinatorial Maps and Linear Cell Complex

-   Added `correct_invalid_attributes`,
    `set_automatic_attributes_management` and
    `are_attributes_automatically_managed` methods in `CombinatorialMap`
    concept. This allows high level operations to not update non void
    attributes during massive calls of these operations, but only after
    the end of their executions.

### 2D Triangulations

-   The class `Constrained_triangulation_plus_2` now can handle
    polylines as constraints.
-   As a consequence a `Constraint_id` has been introduced which
    replaces `pair<Vertex_handle,Vertex_handle>` as identifier of a
    constraint.

### 3D Mesh Generation

-   Add member functions `output_boundary_to_off` and
    `output_facets_in_complex_to_off` in the class
    `CGAL::Mesh_complex_3_in_triangulation_3` to export the boundary of
    a domain or a subdomain.

### 3D Fast Intersection and Distance Computation

-   Add new constructors to `AABB_halfedge_graph_segment_primitive` and
    `AABB_face_graph_triangle_primitive` in order to be able to build
    primitives one by one.

### Spatial Searching

-   Fixed a bug in `CGAL::Splitters.h` sliding midpoint rule, where
    degenerated point sets (e.g.,points on segment) caused the kd-tree
    to get linear.
-   Improved performance of `Orthogonal_k_neighbor_search`. Note that VC
    2013 does not compile `boost::container::deque` of Boost 1\_55 and
    does hence have a workaround which does not have the improvement.
-   **Breaking change:** The concept `OrthogonalDistance` has new
    function overloads for `min_distance_to_rectangle` and
    `max_distance_to_rectangle` with an additional reference parameter
    `std::vector`.
-   **Breaking change:** The order of the points in the iterator range
    `[tree.begin(),tree.end()]` is not the order of insertion of the
    points into the tree. This was not guaranteed before but might have
    been observed and exploited by users.
-   Derived `kd_tree_leaf_node` and `kd_tree_internal_node` from
    `kd_tree_node` to save memory.

### Geometric Object Generators

-   Add a new function `random_convex_hull_in_disc_2` that efficiently
    generates a random polygon as the convex hull of uniform random
    points chosen in a disc.

Release 4.5.2
-------------

Release date: February 2015

### General

-   Fix a bug that prevented the compilation with recent versions of
    Boost (&gt;=1.56) when explicit conversions operators (from C++11)
    are supported. That prevented the compilation with Microsoft Visual
    Studio 2013.

### 3D Convex Hulls

-   Fix a non-robust predicate bug that was showing up when input points
    where lexicographically sorted.

### 3D Mesh Generation

-   Fix a bug in the sliver perturbation optimization method. It could
    create some holes on the surface of the mesh.

Release 4.5.1
-------------

Release date: December 2014

### 3D Mesh Generation

-   Fix a bug in the sliver exudation preservation of boundaries.

Release 4.5
-----------

Release date: October 2014

### Installation

-   Changes in the set of supported platforms:
    -   The Microsoft Windows Visual C++ compiler 2008 (VC9) is no
        longer supported since CGAL-4.5.
-   Since CGAL version 4.0, Eigen was the recommended third-party
    library to use with *Planar Parameterization of Triangulated Surface
    Meshes*, *Surface Reconstruction from Point Sets*, *Approximation of
    Ridges and Umbilics on Triangulated Surface Meshes*, and *Estimation
    of Local Differential Properties of Point-Sampled Surfaces*
    packages. From CGAL version 4.5, Taucs, Blas and Lapack are no
    longer supported.
-   CGAL is now compatible with the new CMake version 3.0.

### Triangulated Surface Mesh Deformation (new package)

-   This package allows to deform a triangulated surface mesh under
    positional constraints of some of its vertices without requiring any
    additional structure other than the surface mesh itself. The methods
    provided implements an as-rigid-as-possible deformation. Note that
    the main class name has changed between the 4.5-beta1 and the 4.5
    releases to better match the CGAL naming conventions (from
    `CGAL::Deform_mesh` to `CGAL::Surface_mesh_deformation`).

### CGAL and the Boost Graph Library (major changes)

-   Cleanup of the `HalfedgeGraph` concept. In particular:
    -   Introduction of the notion of `halfedge_descriptor` in the
        specialization of the class `boost::graph_traits`.
    -   Deprecation of `halfedge_graph_traits`.
    -   A model of `HalfedgeGraph` is considered as an undirected graph.
        Thus any call to `edges()` should be replaced by `halfedges()`
        and `num_edges()` now returns the number of (undirected) edges.
    -   **Breaking change:** `is_border_edge` and `is_border_halfedge`
        properties are removed. The free functions `is_border()` and
        `is_border_edge()` should be used instead.
    -   Renaming of `HalfedgeGraph` specific free functions.
-   Introduction of the `FaceGraph` concept.
-   Adaptation of the package *Triangulated Surface Mesh Simplification*
    and of the class `AABB_halfedge_graph_segment_primitive` from the
    package *3D Fast Intersection and Distance Computation* to the API
    change.
-   Update of the package *Triangulated Surface Mesh Segmentation* and
    of the class `AABB_face_graph_triangle_primitive` from the package
    *3D Fast Intersection and Distance Computation* to accept model of
    the newly introduced concepts.
-   Offer *Euler* operations as free functions for models of the graph
    concepts provided by CGAL.
-   Specialization of `boost::graph_traits` for
    `OpenMesh::PolyMesh_ArrayKernelT` as proof of concept. A
    `OpenMesh::PolyMesh_ArrayKernelT` becomes a model of the
    aforementioned concepts when including
    `CGAL/boost/graph/graph_traits_PolyMesh_ArrayKernelT.h`.

### dD Geometry Kernel

-   A new model `Epick_d` of the `Kernel_d` concept is introduced. It
    provides better performance through arithmetic filtering and
    specializations for fixed dimensions. It may not work with compilers
    as old as gcc-4.2, but was tested with gcc-4.4.

### 3D Convex Hulls

-   Clean up the documentation of the concepts

### 2D Arrangements

-   Fixed a bug in removing an unbounded curve (e.g., a ray) from an
    arrangement induced by unbounded curves.

### 2D Snap Rounding

-   Replaced use of private `kd_tree` with CGAL's official `Kd_tree`
    from `Spatial_searching` package; results in a small performance
    gain. Removed the private `kd_tree` package.

### 3D Triangulations

-   Add an experimental parallel version of the Delaunay triangulation
    and the regular triangulation algorithms, which allows parallel
    insertion and removal of point ranges.
-   Add caching of circumcenters to `Regular_triangulation_cell_base_3`.
    The cache value is computed when `cell->circumcenter()` or
    `rt.dual(cell)` functions are called.

### 3D Periodic Triangulations

-   Add a method to locate point with inexact predicates.

### 3D Mesh Generation

-   Add a new constructor for the class `Labeled_mesh_domain_3` which
    takes an `Iso_cuboid_3`.
-   Add a new labeling function wrapper for meshing multi-domain.
-   The meshing functionality in the Qt demos in `demo/Polyhedron/` and
    `demo/Mesh_3/` can now use the handling of 1d-features, that exists
    in CGAL since version 3.8.
-   Add an experimental parallel version of the 3D mesh refinement and
    mesh optimization methods.

### Point Set Processing and Surface Reconstruction from Point Sets

-   The former demo has been removed and is fully merge in the
    Polyhedron demo.

### Point Set Processing

-   Workaround a bug in dijsktra shortest path of boost 1.54 by shipping
    and using the boost header from the 1.55 release. This header will
    be used only if you are using the version 1.54 of boost.

### Triangulated Surface Mesh Simplification

-   **Breaking change:** Due to the cleanup of the concepts of the
    package *CGAL and the Boost Graph Library*, the named parameter
    `edge_is_border_map` has been removed, and the named parameter
    `edge_is_constrained_map` now expects a property map with an edge
    descriptor as key type (vs. halfedge descriptor before).
-   Add some optimization in the code making the implementation faster
    (depending on the cost and the placement chosen). However, for an
    edge which collapse is not topologically valid, the vector of
    vertices of the link provided by its profile might contains
    duplicates, thus also breaking the orientation guarantee in the
    vector. This must not be a problem for users as the edge is not
    collapsible anyway but if it is a absolute requirement for user
    defined cost/placement, defining the macro
    `CGAL_SMS_EDGE_PROFILE_ALWAYS_NEED_UNIQUE_VERTEX_IN_LINK` will
    restore the former behavior.

### dD Spatial Searching

-   Added methods `reserve(size_t size)` and `size_t       capacity()`
    to class `Kd_tree` to allocate memory to store `size` points and to
    report that number (STL compliance).

### STL Extensions for CGAL

-   Add `Compact_container::operator[]`, allowing a direct access to the
    ith element of a compact container.
-   Add `Concurrent_compact_container`, a compact container which allows
    concurrent insertion and removal.

Release 4.4
-----------

Release date: April 2014

### Installation

-   Additional supported platforms:
    -   The Apple Clang compiler version 5.0 is now supported on
        OS X Mavericks.
    -   The Microsoft Windows Visual C++ compiler 2013 (VC12) is now
        supported.

### Triangulated Surface Mesh Segmentation (new package)

-   This package implements the segmentation of triangulated surface
    meshes based on the Shape Diameter Function (SDF). In addition, it
    also provides functions to generate segmentations based on a user
    defined alternative to the SDF.

### Number Types

-   A new class `CGAL::Mpzf` is introduced on some platforms for exact
    ring operations. It is used to improve the speed of the evaluation
    of predicates in degenerate situations.

### 2D and 3D Geometry Kernel

-   Fix a bug introduced in CGAL 4.3 when computing the intersection of
    two 3D triangles.

### 2D Polygon Partitioning

-   Bug fix to make the partition algorithms working with a Lazy kernel
    such as `Exact_predicates_exact_constructions_kernel`.

### 2D Regularized Boolean Set-Operations

-   Fix two memory leaks in `CGAL::General_polygon_set_2`.

### Combinatorial Maps and Linear Cell Complex

-   `null_dart_handle` is no longer a static data member in the
    `CombinatorialMap` concept. This implies to move the following
    methods of `Dart` concept into `CombinatorialMap` concept:
    `is_free`, `highest_nonfree_dimension`, `opposite` and
    `other_extremity`. We also transform the static methods
    `vertex_attribute` and `point` of `Linear_cell_complex` class into
    non static methods. You can define the CGAL\_CMAP\_DEPRECATED macro
    to keep the old behavior.

### 2D Arrangements

-   Revise the API of **polylines**. In particular, *construction* is
    now done using functors and *iteration* is possible only on the
    segments of a polyline.
-   Fix a bug in the *Landmark* point-location strategy.

### 2D Snap Rounding

-   Fix a memory leak

### 2D Triangulations

-   Add different overloads of the function `insert_constraints` that
    inserts a range of points and segments, or a range of segments.
    These functions uses the spatial sorting in order to speed up the
    time needed for the insertion.

### 3D Alpha Shapes

-   Add member functions in `CGAL::Alpha_shape_3` to give access to the
    alpha status of edges and facets (`get_alpha_status())`.
-   Add another filtration method (`filtration_with_alpha_values()`)
    that reports the alpha value at which each face appears in the
    filtration.

### 3D Mesh Generation

-   Fix the access to functions `number_of_facets` and `number_of_cells`
    in `Mesh_complex_3_in_triangulation_3`.
-   Change the internal API of the sliver perturber, to make possible
    for developers to optimize another criterion than the (default)
    minimal dihedral angle. Developers can also define a new
    perturbation vector (for angles we had gradient of squared
    circumradius, gradient of volume, gradient of minimal dihedral
    angle, and random) which is better suitable to optimize their
    criterion.
-   Improve the use of cache values in `Mesh_cell_base_3` to (re)compute
    circumcenters and sliver criterion values only when needed.

### Triangulated Surface Mesh Simplification

-   Fix a bug in the way edges can be marked as non-removable by adding
    a named-parameter `edge_is_constrained_map` to the function
    `edge_collapse`

### dD Spatial Searching

-   Fix a documentation bug: The property map passed as template
    parameter to the classes `Search_traits_adapter` and
    `Distance_adapter` must be a lvalue property map. To avoid incorrect
    usage, a static assertion has been added in the CGAL code to prevent
    the user from instantiating these classes with an incorrect property
    map type.

### CGAL ipelets

-   Better description of the demo ipelets in the user manual
-   New ipelet for pencils of circles
-   New ipelet for hyperbolic geometry in Poincaré model
-   The generator ipelet now generates point in a selected zone
-   Hilbert sort ipelet implements two policies

Release 4.3
-----------

Release date: October 2013

### The CGAL Manual

-   The documentation of CGAL is now generated with Doxygen.

### 2D Periodic Triangulations (new package)

-   This package allows to build and handle triangulations of point sets
    in the two dimensional flat torus. Triangulations are built
    incrementally and can be modified by insertion or removal of
    vertices. They offer point location facilities. The package provides
    Delaunay triangulations and offers nearest neighbor queries and
    primitives to build the dual Voronoi diagrams.

### API Changes

#### 2D and 3D Geometry Kernel

-   The intersection functions and functors used to return a
    `CGAL::Object` in order to deal with the different possible return
    types. However, depending on the arguments it is possible to reduce
    the possible return types to a small set. For this reason and to
    take advantage of the type safety, we decided to use
    `boost::variant` instead of `CGAL::Object`. The `result_of` protocol
    is now even more useful to determine the return type of the
    intersection functions and functors. The change should be relatively
    transparent to the user thanks to the implicit constructor added to
    `CGAL::Object`. However, it is recommended to upgrade your code. The
    previous behavior can be restored by defining the macro
    `CGAL_INTERSECTION_VERSION` to 1.

#### 2D Arrangements

-   The type of the result of point location queries changed to
    `boost::variant` (from `CGAL::Object`). For convenience, the
    previous behavior can be restored by defining the macro
    `CGAL_ARR_POINT_LOCATION_VERSION` to 1.
-   Introduced an optimization for operations on large and dense
    arrangements.

#### 3D Fast Intersection and Distance Computation

-   Following the intersection API change, `Object_and_primitive_id` has
    been replaced by a template class
    `Intersection_and_primitive_id<Query>` to determine the type
    depending on the query object type.

#### CGAL and Boost Property Maps

-   The `key_type` of the property maps provided by CGAL used to be an
    iterator. In order to be more easily reused, the `key_type` has
    been changed to be the `value_type` of the iterator. The packages
    that have been updated to match these changes are **Point Set
    Processing** and **Surface Reconstruction from Point Sets**.
    However, for most users this change should be transparent if the
    default property maps were used. For convenience, the former
    behavior can be enabled by defining the macro
    `CGAL_USE_PROPERTY_MAPS_API_V1`.

### Algebraic Foundations

-   For convenience, add an overload of `make_rational()` taking a pair
    of numbers.

### 2D and 3D Geometry Kernel

-   A `Iso_rectangle_2` can now be constructed from a `Bbox_2` and an
    `Iso_cuboid_3` from a `Bbox_3`.
-   The implementation of `CGAL::Object` has been updated and now uses
    `boost::shared_ptr` and `boost::any`. This implementation is faster.
-   Add to `Bbox_2` and `Bbox_3` a `+=` operator as well as free
    functions to get the bounding box of a range of geometric objects.

### Combinatorial Maps

-   Two bug fixes: do not use the 2 least significant bits for cell
    attribute without dart support; share the mark when copying a
    CMap\_cell\_iterator.
-   Add a constructor taking a given combinatorial map as argument,
    possibly with different dimension and/or different attributes. This
    allows to transform a combinatorial map.
-   Add operator= and swap method.
-   Add dynamic onmerge/onsplit functions that can be associated
    dynamically to i-attributes and which are automatically called when
    i-cells are split/merged.
-   Add a function allowing to reverse the orientation of a
    combinatorial map, and another one to reverse one connected
    component of a combinatorial map.

### 3D Boolean Operations on Nef Polyhedra

-   Bug-fix in IO when using `Lazy_exact_nt` as number type or
    `Exact_predicates_exact_constructions_kernel` as kernel.

### 2D Triangulations

-   Extend the concept `TriangulationDataStructure_2` to require a more
    general `copy_tds` function that allows a copy between TDS of
    different types. The CGAL model has been updated.
-   Add a way to efficiently insert a range of points with information
    into the 2D constrained Delaunay triangulations.

### 3D Triangulations

-   Extend the concept `TriangulationDataStructure_3` to require a more
    general `copy_tds` function that allows a copy between TDS of
    different types. The CGAL model has been updated.
-   Add an advanced function to set the infinite vertex of the
    triangulation for low level operations
-   Fix a bug in the function inserting a range of points with info when
    the `Fast_location` tag is used

### 2D Segment Delaunay Graph

-   Add functions `insert_points` and `insert_segments` to insert a
    range of points and segments. These functions uses the spatial
    sorting in order to speed up the time needed for the insertion. The
    function
    `insert(Input_iterator first, Input_iterator beyond,       Tag_true)`
    has been updated to dispatch the input when possible to these
    functions.

### 2D Apollonius Graphs

-   Modified insertion algorithm so that the code can handle
    pseudo-circles as well.
-   Updated implementation of the vertex conflict predicate by a faster
    version.

### 3D Mesh Generation

-   Speed-up `Mesh_3` and in particular the global optimizers (Lloyd and
    ODT) by introducing a parameter `do_freeze` to prevent from moving
    vertices which would move of very small displacements.
-   Introduce new data structures and options for speed-up and
    compacity. Note that `Compact_mesh_cell_base_3` and
    `Mesh_vertex_base_3` are now our favored implementations of the
    concepts MeshCellBase\_3 and MeshVertexBase\_3.
-   Introduce a new constructor for `Polyhedral_mesh_domain_3` that
    takes a bounding polyhedron to be meshed along with a polyhedral
    surface entirely included in it. This allows the user to mesh a
    polyhedral domain with internal surface(s) which can be
    non-watertight and even non-manifold.
-   Several documentation bug fixes.
-   Provide the ability to plug in custom cell\_base/vertex\_base
    classes into the Mesh\_triangulation\_3 class.

### Triangulated Surface Mesh Simplification

-   Fix a segmentation fault that was happening when some edges of
    length 0 were in the input mesh.

### 3D Fast Intersection and Distance Computation

-   Following the intersection API change, `Object_and_primitive_id` has
    been replaced by a template class
    `Intersection_and_primitive_id<Query>` to determine the type
    depending on the query object type.
-   Introduce the class `AABB_halfedge_graph_segment_primitive`, which
    replaces the class `AABB_polyhedron_segment_primitive` (which is now
    deprecated). The new class is more general and can be used with any
    model of `HalfedgeGraph`.
-   Introduce the class `AABB_face_graph_triangle_primitive` which
    replaces the class `AABB_polyhedron_triangle_primitive` (which is
    now deprecated).
-   Document the classes `AABB_segment_primitive` and
    `AABB_triangle_primitive` that were already used in some examples.
-   Add a generic primitive class `AABB_primitive` that allows to define
    a primitive type by defining only two property maps.
-   Introduce a new concept of primitive `AABBPrimitiveWithSharedData`.
    It allows to have some data shared between the primitives stored in
    a `AABB_tree`. With this you can, for example have a primitive
    wrapping an integer which refers to the position of a geometric
    object in a `std::vector`. Only one reference to this vector will be
    stored in the traits of the tree. The concept `AABBTraits`, its
    model `AABB_traits` and the class `AABB_tree` have been updated
    accordingly. However, everything is backward compatible.
-   Fix a memory leak in the destructor of the class `AABB-tree`

### STL Extensions for CGAL

-   Add to `Dispatch_output_iterator` and
    `Dispatch_or_drop_output_iterator` an operator to accept and
    dispatch a tuple of values.

### Concurrency in CGAL

-   Add a `FindTBB` CMake module so that one can easily link with TBB to
    write shared-memory parallel code.
-   Introduce two new tags: Sequential\_tag and Parallel\_tag

Release 4.2
-----------

Release date: March 2013

### Installation

-   Additional supported platforms:
    -   The Microsoft Windows Visual C++ compiler 2012 (VC11) is now
        supported.
-   With Microsoft Visual C++ (all supported versions), the compiler
    flags `/bigobj` and `/wd4503` are added by CGAL CMake scripts.
-   This is the last release whose "`UseCGAL.cmake`" file (if using CGAL
    in a CMake build environment) contains the line

          link_libraries(${CGAL_LIBRARIES_DIR} ${CGAL_3RD_PARTY_LIBRARIES_DIRS})

    as this is a deprecated CMake command. The correct way to link with
    CGAL's libraries (as for required 3rd party libraries) is to use
    '`target_link_libraries`' which specifies for each build target
    which libraries should be linked. The following serves as example:

          find_package(CGAL)
          include(${CGAL_USE_FILE})
          add_executable(myexe main.cpp)
          target_link_libraries(myexe ${CGAL_LIBRARIES}
                                      ${CGAL_3RD_PARTY_LIBRARIES})

    We also expect further changes in CGAL's CMake setup (change of
    variable names, consistency of filename and output, removing
    essential libraries, building executables, removal of
    '`${CGAL_3RD_PARTY_LIBRARIES}`').

### 2D Arrangements

-   Enhanced the 2D-arrangements demonstration program and ported it to
    Qt4. The new demonstration program makes use of the CGAL Graphics
    View framework, in which the 2D primitives are individually
    represented as objects in a scene. (The implementations of several
    demos in CGAL already make use of this framework.) This project was
    carried out as part of the 2012 Google Summer of Code program.
-   Fixed a bug in the Walk-Along-A-Line point location strategy for
    arrangements induced by unbounded curves.

### 2D Circular Geometry Kernel

-   Fix the intersection type computed when intersecting two identical
    circles.
-   Forward correctly the result type of the linear kernel functors

### 2D Triangulations

-   Add mechanism to avoid call stack overflow in
    `Delaunay_triangulation_2` and
    `Constrained_Delaunay_triangulation_2`.
-   Add a constructor for `Regular_triangulation_2` and
    `Delaunay_triangulation_2` from a range of points or a range of
    points with info.

### 2D Voronoi Diagram Adaptor

-   Bug-fix: Add ccb() method in face type as documented.

### 3D Minkowski Sum of Polyhedra

-   Fix a memory leak.

### 3D Fast Intersection and Distance Computation

-   Update requirements of the concepts `AABBTraits` and
    `AABBGeomTraits` to match the implementation of the package.

### Generator

-   Addition of the `Combination_enumerator`

### STL Extensions

-   Introduction of `CGAL::cpp11::result_of` as an alias to the tr1
    implementation from boost of the `result_of` mechanism. When all
    compilers supported by CGAL will have a Standard compliant
    implementation of the C++11 `decltype` feature, it will become an
    alias to `std::result_of`.

### Surface Reconstruction from Point Sets

-   Performance improvements and addition of an option to better
    reconstruct undersampled zones. The poisson reconstruction plugin of
    the Polyhedron demo has an option to switch it on.

Release 4.1
-----------

Release date: October 2012

### Installation

-   Additional supported platforms:
    -   The Apple Clang compiler versions 3.1 and 3.2 are now supported
        on Mac OS X.
-   Improved configuration for essential and optional external third
    party software
-   Added more general script to create CMakeLists.txt files:
    `cgal_create_CMakeLists`
-   Availability tests for C++11 features are now performed with the
    help of [Boost.Config](https://www.boost.org/libs/config). A Boost
    version of 1.40.0 or higher is needed to use C++11 features.

### 2D Arrangement

-   Improved the implementation of the incremental randomized
    trapezoidal decomposition point-location strategy. The new
    implementation enables point location in unbounded arrangements. It
    constructs a search structure of guaranteed linear size with
    guaranteed logarithmic query time.

### 2D Convex Hulls and Extreme Points

-   Speed up the preprocessing stage of the Akl-Toussaint implementation
    (used by the free function `convex_hull_2` when forward iterators
    are provided as input).

### Combinatorial Maps

-   Minor bugfix; replace some functors by methods.

### Linear Cell Complex

-   Improve the demo: add a widget showing all the volumes and an
    operation to create a Menger sponge.

### Kernels

-   All Kernel functors now support the result\_of protocol.

### STL\_Extensions for CGAL

-   The namespace `cpp0x` has been renamed `cpp11`. The old name is
    still available for backward compatibility.

Release 4.0.2
-------------

Release date: Jul 2012

This is a bug fix release. It fixes a bug in the `CMakeLists.txt` for
CGAL-4.0.1, that prevented even building the libraries.

Release 4.0.1
-------------

Release date: Jul 2012

This is a bug fix release. Apart various minor fixes in the
documentation, the following has been changed since CGAL-4.0:

### 2D Voronoi Diagram Adaptor (re-added)

-   The package *2D Voronoi Diagram Adaptor* was temporarily removed
    from the CGAL distribution because of license issues. That package
    is now back into CGAL.

### 2D and 3D Geometry Kernel

-   Fix a bug in the `Segment_3-Triangle_3` intersection function in the
    case the segment is collinear with a triangle edge.
-   Fix a bug in the `Projection_traits_.._3` class in the case a
    segment was parallel to the x-axis.

### Algebraic Kernel

-   Avoid the linking error "duplicate symbols" when two compilation
    units using the algebraic kernel are linked.

### 3D Boolean Operations on Nef Polygons Embedded on the Sphere

-   Fix a memory leak due to the usage of an internal mechanism that has
    been replaced by `boost::any`. This also influences the packages 2D
    Boolean Operations on Nef Polygons, 3D Boolean Operations on Nef
    Polyhedra, Convex Decomposition of Polyhedra, and 3D Minkowski Sum
    of Polyhedra.

### 2D Arrangement

-   Fix several memory leaks.

### 2D Mesh Generation

-   Fix a compilation error in the header
    `<CGAL/Mesh_2/Do_not_refine_edges.h>` when g++ version 4.7 is used.

### Surface Mesh Generation and 3D Mesh Generation

-   Fix an important bug in the `CGAL_ImageIO` library, that could lead
    to wrong result when meshing from a 3D image.
-   Fix the compilation of the demo in `demo/Surface_mesher`, when Boost
    version 1.48 or 1.49 is used.

### Surface Mesh Parameterization

-   Fix a memory leak.
-   Fix a compatibility issue with Eigen-3.1 of `Eigen_solver_traits`.
    This fix also affects the usage of that class in the package
    *Surface Reconstruction from Point Sets*.

Release 4.0
-----------

Release date: March 2012

CGAL 4.0 offers the following improvements and new functionality :

### License Changes

The whole CGAL-3.x series was released under a combination of LGPLv2
(for the foundations of CGAL), and QPL (for the high-level packages).
QPL was the former license of the graphical toolkit Qt, but that license
is not supported by any major free software project. Furthermore, the
terms of the LGPLv2 license are ambiguous for a library of C++
templates, like CGAL.

The CGAL project, driven by the CGAL Editorial Board, has decided to
change the license scheme of CGAL. We increased the major number of the
CGAL version to '4' in order to reflect this license change. The
CGAL-4.x series is released under:

-   LGPLv3+ (that is LGPL *"either version 3 of the License, or (at your
    option) any later version"*), for the foundations of CGAL, instead
    of LGPLv2,
-   GPLv3+ for the high-level packages, instead of QPL.

### General

-   On Windows, CGAL libraries are now built by default as shared
    libraries (also called DLL). To run applications that use .dll files
    of CGAL, you must either copy the .dll files into the directory of
    the application, or add the path of the directory that contains
    those .dll files into the PATH environment variable.
-   On Windows, the CMake scripts of CGAL now search for shared version
    of the Boost libraries. You must ensure that the .dll files of Boost
    are found by the dynamic linker. You can, for example, add the path
    to the Boost .dll files to the PATH environment variable.
-   On Windows, CMake version 2.8.6 or higher is now required.
-   Eigen version 3.1 or later is now the recommended third party
    library to use in *Planar Parameterization of Triangulated Surface
    Meshes*, *Surface Reconstruction from Point Sets*, *Approximation of
    Ridges and Umbilics on Triangulated Surface Meshes*, and *Estimation
    of Local Differential Properties of Point-Sampled Surfaces*
    packages. If you use Eigen you no longer need Taucs, Lapack or Blas
    to use those packages (and any other in CGAL).

### Linear Cell Complex (new package)

-   This package implements linear cell complexes, objects in
    d-dimension with linear geometry. The combinatorial part of objects
    is described by a combinatorial map, representing all the cells of
    the object plus the incidence and adjacency relations between cells.
    Geometry is added to combinatorial maps simply by associating a
    point to each vertex of the map. This data structure can be seen as
    the generalization in dD of the `Polyhedron_3`.

### 2D Voronoi Diagram Adaptor (temporarily removed)

-   As the copyright holder of this package has not granted the right to
    switch from QPL to GPL, this package is removed from the
    distribution. Note that it is "only" an adapter, that is the
    functionality of point/segment/disk Voronoi diagram is offered
    through the Delaunay triangulation, segment Delaunay graph, and
    Apollonius graph.

### AABB Tree

-   Document constness of member functions of the `AABB_tree` class.
-   The class `AABB_tree` is now guaranteed to be read-only thread-safe.
    As usual in CGAL, this small overhead introduced for thread-safety
    can be deactivated by defining `CGAL_HAS_NO_THREADS`.

### 2D Alpha Shapes

-   Add an extra template parameter to the class `Alpha_shape_2` that
    allows a certified construction using a traits class with exact
    predicates and inexact constructions.
-   An object of type `Alpha_shape_2` can now be constructed from a
    triangulation.

### 3D Alpha Shapes

-   Add an extra template parameter to the class `Alpha_shape_3` that
    allows a certified construction using a traits class with exact
    predicates and inexact constructions.

### Geometric Object Generators

-   `Random_points_in_iso_box_d` (deprecated since 3.8) has been
    removed. Use `Random_points_in_cube_d` instead.

### Linear and Quadratic Programming Solver

-   Minor bugfix.

### Spatial Searching

-   The const-correctness of this package have been worked out. The
    transition for users should be smooth in general, however adding few
    const in user code might be needed in some cases.
-   The class `Kd_tree` is now guaranteed to be read-only thread-safe.
    As usual in CGAL, this small overhead introduced for thread-safety
    can be deactivated by defining `CGAL_HAS_NO_THREADS`.
-   Bug-fix in `Orthogonal_incremental_neighbor_search` and
    `Incremental_neighbor_search` classes. Several calls to `begin()`
    now allow to make several nearest neighbor search queries
    independently.

### STL Extension

-   `CGAL::copy_n` is now deprecated for `CGAL::cpp0x::copy_n` which
    uses `std::copy_n`, if available on the platform.
-   `CGAL::successor` and `CGAL::predecessor` are now deprecated for
    `CGAL::cpp0x::next` and `CGAL::cpp0x::prev`. These functions use the
    standard versions if available on the platform. Otherwise,
    `boost::next` and `boost::prior` are used.

### Triangulation\_2

-   Fix a thread-safety issue in `Delaunay_triangulation_2` remove
    functions. As usual in CGAL, the small overhead introduced for
    thread-safety can be deactivated by defining `CGAL_HAS_NO_THREADS`.
-   Add extraction operator for the class `Constrained_triangulation_2`
    (and thus to all inheriting classes).

Release 3.9
-----------

Release date: September 2011

CGAL 3.9 offers the following improvements and new functionality :

### General

-   The class `Root_of_2` is now deprecated. It is recommended to use
    the class `Sqrt_extension` instead.
-   The class `Sqrt_extension` is now used everywhere in CGAL where an
    algebraic number of degree 2 is needed. This change has been done in
    the `Root_of_traits` mechanism (indirectly packages 2D Circular
    kernel and 3D Spherical kernel) and the packages 2D Segment Delaunay
    Graphs and 2D Arrangements.
-   Various fixes in the manual.

### Combinatorial Maps (new package)

-   This package provides a new combinatorial data structure allowing to
    describe any orientable subdivided object whatever its dimension. It
    describes all cells of the subdivision and all the incidence and
    adjacency relations between these cells. For example it allows to
    describe a 3D object subdivided in vertices, edges, faces and
    volumes. This data structure can be seen as the generalization in dD
    of the halfedge data structure.

### 3D Convex Hull (major performance improvement)

-   The quickhull implementation of CGAL (`CGAL::convex_hull_3`) has
    been worked out to provide very better performances.
-   The function `CGAL::convex_hull_3` no longer computes the plane
    equations of the facets of the output polyhedron. However an example
    is provided to show how to compute them easily.
-   A global function `convex_hull_3_to_polyhedron_3` is now provided to
    extract the convex hull of a 3D points set from a triangulation of
    these points.

### dD Spatial Searching (major new feature added)

-   A traits-class and distance adapter that together with a point
    property map, allow to make nearest neighbor queries on keys instead
    of points have been added.
-   Few bug fixes in the documentation have revealed some
    inconsistencies that have been corrected. Two traits class concept
    are now documented (`RangeSearchTraits` and `SearchTraits`). Most
    other changes concerns only classes documented as advanced. One
    issue that user can encounter is due to an additional requirement on
    the nested class `Construct_cartesian_const_iterator_d` defined in
    the concept SearchTraits that must provide a nested type
    `result_type`.

### Spatial Sorting (major new feature added)

-   General dimension is now supported.
-   Hilbert sorting admits now two policies: splitting at median or at
    middle (see user manual).
-   Using a property map, sorting on keys instead of points is now
    easier

### dD Kernel

-   The d-dimensional kernel concept and models have been modified to
    additionally provide two new functors `Less_coordinate_d` and
    `Point_dimension_d`.

### 2D Arrangements

-   A new geometry-traits class that handles rational arcs, namely
    `Arr_rational_function_traits_2`, has been introduced. It replaced
    an old traits class, which handled the same family of curves, but it
    was less efficient. The new traits exploits CGAL algebraic kernels
    and polynomials, which were not available at the time the old traits
    class was developed.
-   A new geometry traits concept called
    `ArrangementOpenBoundaryTraits_2` has been introduced. A model of
    this concept supports curves that approach the open boundary of an
    iso-rectangular area called parameter space, which can be unbounded
    or bounded. The general code of the package, however, supports only
    the unbounded parameter space. We intend to enhance the general code
    to support also bounded parameter spaces in a future release.
-   The deprecated member function `is_at_infinity()` of
    `Arrangement_2::Vertex` has been removed. It has been previously
    replaced new function `is_at_open_boundary()`.
-   The tags in the geometry traits that indicate the type of boundary
    of the embedding surface were replaced by the following new tags:

                 Left_side_category
                 Bottom_side_category
                 Top_side_category
                 Right_side_category

    It is still possible not to indicate the tags at all. Default values
    are assumed. This however will produce warning messages, and should
    be avoided.

Release 3.8
-----------

Release date: April 2011

CGAL 3.8 offers the following improvements and new functionality :

### General

-   Boost version 1.39 at least is now required.
-   Initial support for the LLVM Clang compiler (prereleases of version
    2.9).
-   Full support for the options -strict-ansi of the Intel Compiler 11,
    and -ansi of the GNU g++ compiler.
-   Adding a concept of ranges. In the following releases, it will be
    the way to provide a set of objects (vs. a couple of iterators).
-   Fix a memory leak in CORE polynomials.
-   Various fixes in the manual.

### 3D Mesh Generation (major new feature added)

-   Adding the possibility to handle sharp features: the 3D Mesh
    generation package now offers the possibility to get in the final
    mesh an accurate representation of 1-dimensional sharp features
    present in the description of the input domain.

### 2D Triangulations (major new feature added)

-   Add a way to efficiently insert a range of points with information
    into a 2D Delaunay and regular triangulation.
-   Add member function mirror\_edge taking an edge as parameter.
-   Fix an infinite loop in constrained triangulation.

### 3D Triangulations (major new feature added)

-   Add a way to efficiently insert a range of points with information
    into a 3D Delaunay and regular triangulation.
-   Add a member function to remove a cluster of points from a Delaunay
    or regular triangulation.
-   function vertices\_in\_conflict is renamed
    vertices\_on\_conflict\_zone\_boundary for Delaunay and regular
    triangulation. Function vertices\_inside\_conflict\_zone is added to
    regular triangulation.
-   Structural filtering is now internally used in locate function of
    Delaunay and regular triangulation. It improves average construction
    time by 20%.
-   Added demo.

### 3D Alpha Shapes (major new feature added)

-   The new class Fixed\_alpha\_shape\_3 provides a robust and faster
    way to compute one alpha shape (with a fixed value of alpha).

### AABB tree

-   Adding the possibility to iteratively add primitives to an existing
    tree and to build it only when no further insertion is needed.

### 2D and 3D Kernel

-   Better handling of 2D points with elevation (3D points projected
    onto trivial planes). More general traits classes
    (Projection\_traits\_xy\_3,
    Projection\_traits\_yz\_3,Projection\_traits\_yz\_3) are provided to
    work with triangulations, algorithms on polygons, alpha-shapes,
    convex hull algorithm... Usage of former equivalent traits classes
    in different packages is now deprecated.
-   Exact\_predicates\_exact\_constructions\_kernel now better use the
    static filters which leads to performance improvements.
-   Add an overload for the global function angle, taking three 3D
    points.
-   In the 2D and 3D kernel concept, the constant Boolean
    Has\_filtered\_predicates is now deprecated. It is now required to
    use Has\_filtered\_predicates\_tag (being either Tag\_true or
    Tag\_false).
-   Compare\_distance\_2 and Compare\_distance\_3 provide additional
    operators for 3 and 4 elements.
-   Add intersection test and intersection computation capabilities
    between an object of type Ray\_3 and either an object of type
    Line\_3, Segment\_3 or Ray\_3.
-   Improve intersection test performance between an object of type
    Bbox\_3 and an object of type Plane\_3 or Triangle\_3 by avoiding
    arithmetic filter failures.

### 2D Envelope

-   Env\_default\_diagram\_1 is deprecated, Envelope\_diagram\_1 should
    be used instead.

### 3D Envelope

-   A new demo program called `L1_Voronoi_diagram_2` has been
    introduced. It demonstrates how 2D Voronoi diagrams of points under
    the L1 metric are constructed using lower envelopes.

### dD Kernel

-   Add functor Compute\_coordinate\_d to Kernel\_d concept.

### Geometric Object Generators

-   CGAL::Random uses boost::rand48 instead of std::rand.
-   Adding to CGAL::Random a way to generate random integers.
-   Adding generators for dD points.

### Algebraic Foundations

-   Algebraic\_structure\_traits now provides an Inverse functor for
    Fields. There is also a new global function inverse.

### Bounding Volumes

-   dD Min sphere of spheres has a new traits class for the min sphere
    of points.

### Triangulated Surface Mesh Simplification

-   The priority queue internally used to prioritize edge
    simplifications is no longer a relaxed heap but a binomial heap.
    This fix guarantees that all edges satisfying a simplification
    criteria are removed (if possible).

### 3D Boolean Operations on Nef Polyhedra

-   Allow construction of a 3D nef polyhedron from a 3D polyhedron with
    normals.

### 2D Arrangements

-   Fix a bug in the method insert\_at\_vertices of the Arrangement\_2
    class.
-   Fix several bugs in the traits class Arr\_Bezier\_curve\_traits\_2
    for arrangement of Bezier curves.

### 2D Minkowski Sums

-   A bug in the convolution method was fixed.

Release 3.7
-----------

Release date: October 2010

CGAL 3.7 offers the following improvements and new functionality :

### General

-   The configuration of CGAL libraries now requires CMake&gt;=2.6.
-   Changes in the set of supported platforms:
    -   GNU g++ 4.5 supported (with or without the compilation option
        -std=c++0x).
    -   Initial support for the option -strict-ansi of the Intel
        Compiler 11. The CGAL libraries compile with that option, and
        most CGAL headers have been fixed. The packages "3D Boolean
        Operations on Nef Polyhedra" (Nef\_3), "Convex Decomposition of
        Polyhedra" (Convex\_decomposition\_3), and "3D Minkowski Sum of
        Polyhedra" (Minkowski\_sum\_3) are known to still fail to
        compile with that compiler flag.
    -   The Microsoft Windows Visual C++ compiler 2010 (VC10), that was
        experimentally supported by CGAL-3.6.1, is now fully supported.
        Note that CMake&gt;=2.8.2 is required for that support.
    -   The Microsoft Windows Visual C++ compiler 2005 (VC8) is no
        longer supported by the CGAL project since CGAL-3.7.
    -   With Microsoft Windows Visual C++ (VC9 and VC10), the optional
        dependencies Gmp, Mpfr, Blas, Lapack, Taucs no longer use
        Boost-style name mangling. Only one variant is now provided by
        the CGAL Windows installer (release, with dynamic runtime).
-   Some demos now require a version of Qt4 &gt;= 4.3.
-   CGAL\_PDB is no longer provided with CGAL. An alternative solution
    for people interested in reading PDB files is to use ESBTL
    (https://esbtl.sourceforge.net/).
-   Fix issues of the CGAL wrappers around the CORE library, on 64 bits
    platforms.

### Arithmetic and Algebra

-   New models Algebraic\_kernel\_d\_1 and Algebraic\_kernel\_d\_2 for
    the corresponding concepts. They provide generic support for various
    coefficient types

### Arrangements

-   A new model Arr\_algebraic\_segment\_traits\_2 of
    ArrangementTraits\_2 that supports algebraic curves of arbitrary
    degree in the plane

### 2D Triangulations

-   The Delaunay and regular 2D triangulations now use a symbolic
    perturbation to choose a particular triangulation in co-circular
    cases.
-   The return type of the template member function insert(It beg, It
    end), taking an iterator range of points, has been changed from int
    to std::ptrdiff\_t.
-   Classes Triangulation\_euclidean\_traits\_xy\_3,
    Triangulation\_euclidean\_traits\_yz\_3 and
    Triangulation\_euclidean\_traits\_xz\_3 are now model of the concept
    ConstrainedTriangulationTraits\_2. They can be used with and without
    intersection of constraints.
-   2D Delaunay and basic triangulations now provide vertex relocation
    by the mean of these two new methods: move and
    move\_if\_no\_collision. The methods are also available for the
    hierarchy (Triangulation\_hierarchy\_2).

### 3D Triangulations

-   The return type of the template member function insert(It beg, It
    end), taking an iterator range of points, has been changed from int
    to std::ptrdiff\_t.
-   3D Delaunay triangulations now provide vertex relocation by the mean
    of these two new methods: move and move\_if\_no\_collision. This
    works in both Compact\_policy and Fast\_policy.

### 2D and 3D Alpha Shapes

-   The type int in the API has been changed to std::size\_t so that
    CGAL can deal with large data sets (64 bit addresses).

### 2D Mesh Generation

-   The execution of the 2D mesh generator is now deterministic (same at
    each run).

### 3D Mesh Generation

-   The efficiency of the 3D mesh generator has been improved (the
    number of calls to the oracle per inserted vertex has globally
    decrease). This is achieved through a slight change of the mesh
    generator strategy which implies that a surface component that is
    not detected at the surface mesher level will never be discovered by
    chance, owing to the refinement of some tetrahedra, as it could
    happen before. Please note that defining the macro
    CGAL\_MESH\_3\_USE\_OLD\_SURFACE\_RESTRICTED\_DELAUNAY\_UPDATE
    switches back to the old behavior.
-   A demo program is now available.

### Surface Reconstruction from Point Sets

-   Improved performance and minor bug fix.

### 2D Range and Neighbor Search

-   The type int in the API has been changed to std::size\_t so that
    CGAL can deal with large data sets (64 bit addresses).

Release 3.6.1
-------------

Release date: June 2010

This is a bug fix release. The following has been changed since
CGAL-3.6:

### General

-   Fix compilation errors with recent Boost versions (since 1.40).
-   Initial support for the Microsoft Visual C++ compiler 10.0 (MSVC
    2010). For that support, CMake&gt;=2.8.2 is required. Note also that
    the compiler option "/bigobj" is necessary to compile some CGAL
    programs with MSVC 2010.

### Polynomial

-   Fix compilation errors with the Microsoft Visual C++ compiler and
    the Intel C++ compiler.

### Polyhedron

-   Fix a compilation errors in demo/Polyhedron/:
-   issue with the location of qglobal.h of Qt4 on MacOS X,
-   missing texture.cpp, if TAUCS is used,
-   Fix the location of built plugins of demo/Polyhedron/, when CGAL is
    configured with WITH\_demos=ON

### 3D Periodic Triangulations

-   Fixed bug in the triangulation hierarchy for periodic
    triangulations.

### 2D Mesh Generation

-   Fix a bug that lead to precondition violation.
-   Improve the user manual about the member function is\_in\_domain()
    of the Face type.
-   The 2D meshing process is now deterministic (sorting of bad faces no
    longer relies on pointers comparisons).

### 3D Mesh Generation

-   Fix a linking errors (duplicate symbols) when
    `<CGAL/refine_mesh_3.h>` is included in different compilation units.

### Spatial Searching

-   Fix a bug in `<CGAL/Orthogonal_k_neighbor_search.h>` when several
    nearest neighbors are at the same distance from the query point.

### IO Streams

-   Fix a bug in `<CGAL/IO/VRML_2_ostream.h>` that generated VRML 2
    files with an invalid syntax for IndexedFaceSet nodes.

### Triangulation\_2

-   Add missing Compare\_distance\_2 functor in trait classes
    Triangulation\_euclidean\_traits\_xy\_3
    Triangulation\_euclidean\_traits\_yz\_3 and
    Triangulation\_euclidean\_traits\_xz\_3. This was preventing calling
    member function nearest\_vertex of Delaunay\_triangulation\_2
    instantiated with one of these traits.

Release 3.6
-----------

Release date: March 2010

CGAL 3.6 offers the following improvements and new functionality :

### General

-   Boost version 1.34.1 at least is now required.

### Arithmetic and Algebra

#### Algebraic Kernel (new package)

-   This new package is targeted to provide black-box implementations of
    state-of-the-art algorithms to determine, compare and approximate
    real roots of univariate polynomials and bivariate polynomial
    systems. It includes models of the univariate algebraic kernel
    concept, based on the library RS.

#### Number Types

-   Two new arbitrary fixed-precision floating-point number types have
    been added: the scalar type Gmpfr and the interval type Gmpfi, based
    on the MPFR and MPFI libraries respectively.

### Geometry Kernels

#### 2D and 3D Geometry Kernel

-   Add new do\_intersect() and intersection() overloads:
    -   do\_intersect(Bbox\_3, Bbox\_3/Line\_3/Ray\_3/Segment\_3)
    -   intersection(Triangle\_3, Line\_3/Ray\_3/Segment\_3)

### Polygons

#### 2D Regularized Boolean Set-Operations

-   Fixed General\_polygon\_set\_2::arrangement() to return the proper
    type of object.

### Arrangement

#### 2D Arrangements

-   Fixed passing a (const) traits object to the constructor of
    Arrangement\_2.
-   Introduced Arrangement\_2::fictitious\_face(), which returns the
    fictitious face in case of an unbounded arrangement.
-   Fixed a bug in Bezier-curve handling.
-   Added (back) iterator, number\_of\_holes(), holes\_begin(), and
    holes\_end() to the default DCEL for backward compatibility.
-   Added (simple) versions of the free overlay() function. It employs
    the default overlay-traits, which practically does nothing.

### Polyhedron

-   Fix a compilation errors in demo/Polyhedron/:
    -   issue with the location of qglobal.h of Qt4 on MacOS X,
    -   missing texture.cpp, if TAUCS is used,
-   Fix the location of built plugins of demo/Polyhedron/, when CGAL is
    configured with WITH\_demos=ON
-   Fix a bug in test\_facet function of the incremental builder: the
    function did not test if while a new facet makes a vertex manifold,
    no other facet incident to that vertex breaks the manifold property.

### Triangulations and Delaunay Triangulations

#### 2D/3D Regular Triangulations

-   Weighted\_point now has a constructor from Cartesian coordinates.

#### 3D Triangulations

-   Regular\_triangulation\_3 : semi-static floating-point filters are
    now used in its predicates, which can speed up its construction by a
    factor of about 3 when
    Exact\_predicates\_inexact\_constructions\_kernel is used.
-   The class Regular\_triangulation\_filtered\_traits\_3 is deprecated,
    the class Regular\_triangulation\_euclidean\_traits\_3 must be used
    instead. The predicates of that traits will be filtered if the
    kernel given as template parameter of that traits is itself a
    filtered kernel.
-   Triangulation\_hierarchy\_3 is now deprecated, and replaced by a
    simpler CGAL::Fast\_location policy template parameter of
    Delaunay\_triangulation\_3.
-   The old version of remove() (enabled with
    CGAL\_DELAUNAY\_3\_OLD\_REMOVE) has been deleted.

#### 3D Periodic Triangulations

-   New demo: 3D periodic Lloyd algorithm.
-   New functionality for Voronoi diagrams: dual of an edge and of a
    vertex, volume and centroid of the dual of a vertex.
-   The package can now be used with the 3D Alpha Shapes package to
    compute periodic alpha shapes.

#### 3D Alpha shapes

-   The class Weighted\_alpha\_shape\_euclidean\_traits\_3 is
    deprecated, the class Regular\_triangulation\_euclidean\_traits\_3
    must be used instead.
-   The package can now be used together with the 3D Periodic
    Triangulation package to compute periodic alpha shapes.

#### 2D/3D Triangulations, 2D Segment Delaunay Graph, 2D Apollonius Graph, and 3D Periodic Triangulations

-   The constructor and insert function taking ranges now produce
    structures whose iterator orders is now deterministic (same at each
    run).

### Mesh Generation

#### 2D Mesh Generation

-   The 2D mesh generator can now be used with a constrained Delaunay
    triangulation with constraints hierarchy
    (Constrained\_triangulation\_plus\_2).
-   In some cases (refinement of a constrained edge that is on the
    convex hull), the 2D mesh generator from CGAL-3.4 and CGAL-3.5 could
    create invalid triangulations. This bug is now fixed.

#### 3D Mesh Generation

-   The mesh generator has been enriched with an optimization phase to
    provide 3D meshes with well shaped tetrahedra (and in particular no
    slivers). The optimization phase involves four different
    optimization processes: two global optimization processes (ODT and
    Lloyd), a perturber and an exuder. Each of these processes can be
    activated or not, and tuned to the users needs and to available
    computer resources.

### Support library

#### CGAL ipelets

-   Add support for version 7 of Ipe.

Release 3.5.1
-------------

Release date: December 2009

This is a bug fix release.

### Documentation

-   Fixes in the documentation (the online documentation of CGAL-3.5 is
    now based on CGAL-3.5.1).
-   Fixes to the bibliographic references.

### Windows installer

-   The Windows installer of CGAL-3.5.1 fixes an issue with downloading
    of precompiled binaries of the external library TAUCS.

### Bug fixes in the following CGAL packages

#### AABB tree

-   Fix a linker issue in do\_intersect(Bbox\_3,Bbox\_3).
-   Fix compilation issue in do\_intersect(Bbox\_3,Ray\_3) when using
    the parameters in this order.

#### 3D Mesh Generation

-   Fix a bug in initial points construction of a polyhedral surface.

Release 3.5
-----------

Release date: October 2009

CGAL releases will now be published about every six months. As a
transition release, CGAL-3.5 has been developed during 9 months from the
release CGAL-3.4.

Version 3.5 differs from version 3.4 in the platforms that are supported
and in functionality. There have also been a number of bug fixes for
this release.

### General

-   Additional supported platforms:
    -   GNU g++ 4.4 supported.
    -   Intel Compiler 11 supported on Linux
-   Fixed ABI incompatibilities when mixing CGAL and Boost Program
    Options on Windows/Visual C++ (the compilation flag
    -D\_SECURE\_SCL=0 is not longer use in Debug mode).

### Geometry Kernels

#### 3D Spherical Geometry Kernel

-   Add functionalities to manipulate circles, circular arcs and points
    that belong to the same sphere.

### Polygons

#### 2D Regularized Boolean Set-Operations

-   The polygon validation operations were enhanced and their interface
    was improved. They are now offered as free functions and applied
    properly.

#### 2D Straight Skeleton and Polygon Offsetting

-   Updated the manual to document the new partial skeletons feature
    (already in the code since 3.4)

### Arrangements

#### 2D Arrangements

-   The member function is\_at\_infinity() of Arrangement\_2::Vertex was
    replaced by the new function is\_at\_open\_boundary(). The former is
    deprecated. While still supported in version 3.5, It will not be
    supported in future releases. The member functions
    boundary\_type\_in\_x() and boundary\_type\_in\_y() were permanently
    replaced by the functions parameter\_space\_in\_x() and
    parameter\_space\_in\_y(), respectively. The 2 new functions return
    an enumeration of a new type, namely Arr\_parameter\_space.
-   The tags in the geometry traits that indicate the type of boundary
    of the embedding surface were replaced by the following new tags:
    Arr\_left\_side\_tag Arr\_bottom\_side\_tag Arr\_top\_side\_tag
    Arr\_right\_side\_tag In addition, the code was change, and now it
    is possible not to indicate the tags at all. Default values are
    assumed. This however will produce warning messages, and should be
    avoided.
-   All operations of the geometry traits-class were made 'const'. This
    change was reflected in the code of this package and all other
    packages that are based on it. Traits classes that maintain state,
    should declare the data members that store the state as mutable.

#### Envelopes of Surfaces in 3D

-   A few bugs in the code that computes envelopes were fixed, in
    particular in the code that computes the envelopes of planes.

### Triangulations and Delaunay Triangulations

#### 3D Periodic Triangulations (new package)

-   This package allows to build and handle triangulations of point sets
    in the three dimensional flat torus. Triangulations are built
    incrementally and can be modified by insertion or removal of
    vertices. They offer point location facilities.

### Mesh Generation

#### Surface Reconstruction from Point Sets (new package)

-   This CGAL package implements an implicit surface reconstruction
    method: Poisson Surface Reconstruction. The input is an unorganized
    point set with oriented normals.

#### 3D Mesh Generation (new package)

-   This package generates 3 dimensional meshes. It computes isotropic
    simplicial meshes for domains or multidomains provided that a domain
    descriptor, able to answer queries from a few different types on the
    domain, is given. In the current version, Mesh\_3 generate meshes
    for domain described through implicit functional, 3D images or
    polyhedral boundaries. The output is a 3D mesh of the domain volume
    and conformal surface meshes for all the boundary and subdividing
    surfaces.

### Geometry Processing

#### Triangulated Surface Mesh Simplification

-   BREAKING API change in the passing of the visitor object.
-   Fixed a bug in the link\_condition test
-   Added a geometric test to avoid folding of facets
-   Fixed a bug in the handling of overflow in the LindstromTurk
    computations
-   Updated the manual to account for the new visitor interface

#### Point Set Processing (new package)

-   This packages implements a set of algorithms for analysis,
    processing, and normal estimation and orientation of point sets.

### Spatial Searching and Sorting

#### AABB tree (new package)

-   This package implements a hierarchy of axis-aligned bounding boxes
    (a AABB tree) for efficient intersection and distance computations
    between 3D queries and sets of input 3D geometric objects.

### Support Library

#### CGAL\_ipelets (new package):

-   Object that eases the writing of Ipe's plugins that use CGAL.
    Plugins for CGAL main 2D algorithm are provided as demo.

Release 3.4
-----------

Release date: January 2009

Version 3.4 differs from version 3.3.1 in the platforms that are
supported and in functionality. There have also been a number of bug
fixes for this release.

### General

-   GNU g++ 4.3 supported. Support for g++ 3.3 is dropped.
-   Visual 9 supported. Support for Visual 7 is dropped.
-   Boost version 1.33 at least is now required.
-   CGAL now depends on Boost.Threads, which implies to link against a
    compiled part of Boost.
-   The new macro CGAL\_NO\_DEPRECATED\_CODE can be defined to disable
    deprecated code, helping users discover if they rely on code that
    may be removed in subsequent releases.
-   Assertion behavior: It is not possible anymore to set the CONTINUE
    mode for assertion failures. Functions that allow to change the
    assertion behavior are now declared in
    `<CGAL/assertions_behaviour.h>`.
-   Qt3 based demos are still there but the documentation has been
    removed as the CGAL::Qt\_Widget will be deprecated.
-   Qt4 based demos use the Qt GraphicsView framework and the
    libQGLViewer.

### Installation

-   install\_cgal has been replaced by CMake.

### Polynomial (new package)

-   This package introduces a concept Polynomial\_d, a concept for
    multivariate polynomials in d variables.

### Modular Arithmetic (new package)

-   This package provides arithmetic over finite fields.

### Number Types

-   the counter Interval\_nt::number\_of\_failures() has been removed,
    replaced by a profiling counter enabled with CGAL\_PROFILE.
-   Fix of a bug in CORE/Expr.h; as a consequence, the arrangement demo
    works properly when handling arrangements of conics, for example,
    when defining an arc with 5 points.

### 3D Spherical Geometry Kernel (new package)

-   This package is an extension of the linear CGAL Kernel. It offers
    functionalities on spheres, circles, circular arcs and line segments
    in the 3D space.

### Linear Kernel

-   We recommend that you use the object\_cast() function instead of
    assign() to extract an object from a CGAL::Object, for efficiency
    reasons.
-   The Kernel archetypes provided by the 2D/3D linear kernel have been
    removed.
-   The deprecated linear kernel functors Construct\_supporting\_line\_2
    and Construct\_supporting\_line\_3 have been removed.
-   Ambiant\_dimension and Feature\_dimenison have been added to
    retrieve the potentially compile-time dimension of a space or of an
    object.
-   barycenter() functions have been added.
-   The geometric object Circle\_3 as well as predicates and
    constructions have been added to the kernel
-   The missing intersection/do\_intersect between Line\_3 and Line\_3
    has been added as well.

### 3D Triangulations

-   Removed the deprecated functions Cell:mirror\_index() and
    Cell::mirror\_vertex().
-   Derecursification of two functions that in some cases lead to stack
    overflows

### 3D Nef Polyhedron

-   n-ary union/intersection
-   intersection with halfspace under standard kernel
-   constructor for polylines

### CGAL and the Qt4 GraphicsView (new package)

-   2D CGAL Kernel objects and many data structures have can be rendered
    in a QGraphicsView

### STL Extensions:

-   The functor adaptors for argument binding and composition (bind\_\*,
    compose, compose\_shared, swap\_\*, negate, along with the helper
    functions set\_arity\_\* and Arity class and Arity\_tag typedefs)
    which were provided by `<CGAL/functional.h>` have been removed.
    Please use the better boost::bind mechanism instead. The concept
    AdaptableFunctor has been changed accordingly such that only a
    nested result\_type is required.
-   The accessory classes Twotuple, Threetuple, Fourtuple and Sixtuple
    are also deprecated (use CGAL::array instead).
-   CGAL::Triple and CGAL::Quadruple are in the process of being
    replaced by boost::tuple. As a first step, we strongly recommend
    that you replace the direct access to the data members (.first,
    .second, .third, .fourth), by the get&lt;i&gt;() member function;
    and replace the make\_triple and make\_quadruple maker functions by
    make\_tuple.
    This way, in a further release, we will be able to switch to
    boost::tuple more easily.
-   The class CGAL::Uncertain&lt;&gt; has been documented. It is
    typically used to report uncertain results for predicates using
    interval arithmetic, and other filtering techniques.

### 2D Arrangements

-   Changed the name of the arrangement package from Arrangement\_2 to
    Arrangement\_on\_surface\_2 to reflect the potential capabilities of
    the package to construct and maintain arrangements induced by curves
    embedded on two dimensional surfaces in three space. Most of these
    capabilities will become available only in future releases though.
-   Enhanced the geometry traits concept to handle arrangements embedded
    on surfaces. Each geometry-traits class must now define the
    'Boundary\_category' tag.
-   Fixed a bug in Arr\_polyline\_traits\_2.h, where the operator that
    compares two curves failed to evaluate the correct result (true)
    when the curves are different, but their graphs are identical.
-   Permanently removed IO/Arr\_postscript\_file\_stream.h and
    IO/Polyline\_2\_postscript\_file\_stream.h, as they depend on
    obsolete features and LEDA.
-   Fixed several bugs in the arrangement demo and enhanced it. e.g.,
    fixed background color change, allowed vertex coloring , enabled
    "smart" color selection, etc.
-   Enhanced the arrangement demo with new features, such as allowing
    the abortion of the merge function (de-select), updated the how-to
    description, etc.
-   Replace the functions CGAL::insert\_curve(), CGAL::insert\_curves(),
    CGAL::insert\_x\_monotone\_curve(), and
    CGAL::insert\_x\_monotone\_curves() with a single overloaded
    function CGAL::insert(). The former 4 functions are now deprecated,
    and may no longer be supported in future releases.

### Envelopes of Surfaces in 3D

-   Fixed a bug in the computation of the envelope of unbounded planes
    caused by multiple removals of vertices at infinity.

### 2D Regularized Boolean Set-Operations

-   Fixed a bug in connect\_holes() that caused failures when connecting
    holes touching the outer boundary.
-   Fixed the concept GeneralPolygonSetTraits\_2. Introduced two new
    concepts GpsTraitsGeneralPolygon\_2 and
    GpsTraitsGeneralPolygonWithHoles\_2. Fixed the definition of the two
    nested required types Polygon\_2 and Polygon\_with\_holes\_2 of the
    GeneralPolygonSetTraits\_2 concept. They must model now the two new
    concepts above.
-   Added a default template parameter to 'General\_polygon\_set\_2' to
    allow users to pass their specialized DCEL used to instantiate the
    underlying arrangement.
-   Enhanced the BOP demo to use multiple windows.

### 2D Minkowski Sums

-   Fixed a few bugs in the approximate offset function, making it
    robust to highly degenerate inputs.
-   Fixed a bug in the exact Minkowski sum computation when processing
    degenerate inputs that induce overlapping of contiguous segments in
    the convolution cycles.
-   Optimized the approximate offset function (reduced time consumption
    up to a factor of 2 in some cases).
-   Added functionality to compute the offset (or to approximate the
    offset) of a Polygon\_with\_holes\_2 (and not just of a Polygon\_2).
-   Added the functionality to compute (or to approximate) the inner
    offset of a polygon.

Release 3.3.1
-------------

Release date: August 2007

This is a bug fix release.

### General

-   Intel C++ 9 was wrongly recognized as unsupported by install\_cgal.
-   Added autolink (for Visual C++) for the CGALImageIO and CGALPDB
    libraries.
-   Fixed bug in Memory\_sizer when using more than 4GB of memory
    (64bit).

### Number Types

-   Fixed bug in FPU rounding mode macros (affected only the alpha
    architecture).
-   Fixed bug in MP\_Float constructor from double for some particular
    values.
-   Fixed bug in to\_double(Lazy\_exact\_nt) sometimes returning NaN.

### Kernel

-   Fixed forgotten derivation in Circular\_kernel\_2::Has\_on\_2
-   Added some missing functions in Bbox\_3 compared to Bbox\_2.

### Skin Surface Meshing

-   The new Skin Surface Meshing package had been forgotten in the list
    of changes and the release announcement of CGAL 3.3: This package
    allows to build a triangular mesh of a skin surface. Skin surfaces
    are used for modeling large molecules in biological computing.

### Arrangements

-   Fixed a bug in the Arrangement\_2 package in dual arrangement
    representation for Boost graphs when reporting all halfedges of a
    face.
-   Fixed a bug in the Arrangement sweep-line when using a specific
    polyline configuration.
-   Fixed bug in Arrangement\_2 in walk along a line point location for
    unbounded curves.
-   Fixed bug in aggregated insertion to Arrangement\_2.
-   Fixed bug in Arrangement\_2 class when inserting an unbounded curve
    from an existing vertex.
-   Fixed bug when dealing with a degenerate conic arc in
    Arr\_conic\_traits\_2 of the Arrangement package, meaning a line
    segment which is part of a degenerate parabola/hyperbola.
-   Fixed bug in the Bezier traits-class: properly handle line segments.
    properly handle comparison near a vertical tangency.

### 2D Polygon

-   Fixed bug in degenerate case of Polygon\_2::is\_convex() for equal
    points.

### 2D Triangulations

-   Fixed bug in Regular\_triangulation\_2.

### 3D Triangulations

-   Added a circumcenter() function in the default Cell type parameter
    Triangulation\_ds\_cell\_base\_3, so that the .dual() member
    function of Delaunay still work as before, without requiring the
    explicit use of Triangulation\_cell\_base.
-   Added missing operator-&gt;() to Facet\_circulator.

### Interpolation

-   Fixed bug in Interpolation 3D about the normalization coefficient
    initialization.

### 3D Boolean Operations on Nef Polyhedra

-   Fixed bug in construction of Nef\_polyhedron\_3 from off-file. Now,
    always the inner volume is selected.
-   Fixed bug in conversion from Nef\_polyhedron\_3 to Polyhedron\_3.
    Polyhedron\_3 was not cleared at the beginning.
-   Fixed bug in Nef\_polyhedron\_3 in update of indexes for
    construction of external structure.

### Third Party Libraries Shipped with CGAL

-   TAUCS supports now 64 bits platforms.
-   CAUTION: Since version 3.3.1, CGAL is no longer compatible with the
    official release of TAUCS (currently 2.2). Make sure to use the
    version modified by the CGAL project and available from the download
    section of https://www.cgal.org.

Release 3.3
-----------

Release date: May 2007

Version 3.3 differs from version 3.2.1 in the platforms that are
supported and in functionality. There have also been a number of bug
fixes for this release.

Additional supported platforms

-   GNU g++ 4.1 and 4.2
-   Intel C++ compiler 9
-   Microsoft Visual C++ compiler 8.0

The following platforms are no longer supported:

-   Intel 8

CGAL now supports Visual C++ "Checked iterators" as well as the debug
mode of GNU g++'s STL (-D\_GLIBCXX\_DEBUG).

CGAL now works around the preprocessor macros 'min' and 'max' defined in
`<windows.h>` which were clashing with min/max functions.

### Installation

-   On Windows the libraries built in Developer Studio now have names
    which encode the compiler version, the runtime and whether it was
    built in release or debug mode. The libraries to link against are
    chosen with linker pragmas in header files.
-   On all platforms but Windows shared and static versions of the
    libraries are generated

### Manuals

-   The Package Overview page now also hosts the precompiled demos.

### Algebraic Foundations

-   Algebraic Foundations (new package)
    This package defines what algebra means for CGAL, in terms of
    concepts, classes and functions. The main features are: (i) explicit
    concepts for interoperability of types (ii) separation between
    algebraic types (not necessarily embeddable into the reals), and
    number types (embeddable into the reals).
-   Number Types
    Fixed\_precision\_nt and Filtered\_exact number types have been
    removed.

### Kernels

-   2D Circular Kernel
    Efficiency improved through geometric filtering of predicates,
    introduced with the filtered kernel
    Filtered\_bbox\_circular\_kernel\_2&lt;.&gt;, and also chosen for
    the predefined kernel Exact\_circular\_kernel\_2.
-   Linear Kernel
    Exact\_predicates\_exact\_constructions\_kernel memory and run-time
    improvements through usage of lazy geometric constructions instead
    of lazy arithmetic.

### Data Structures and Algorithms

-   Surface Mesh Simplification (new package)
    This package provides a mesh simplification framework using edge
    collapse operations, and provides the Turk/Lindstrom simplification
    algorithm.
-   Skin Surface Meshing (new package)
    This package allows to build a triangular mesh of a skin surface.
    Skin surfaces are used for modeling large molecules in biological
    computing. The surface is defined by a set of balls, representing
    the atoms of the molecule, and a shrink factor that determines the
    size of the smooth patches gluing the balls together.
-   Estimation of Local Differential Properties (new package)
    This package allows to compute local differential quantities of a
    surface from a point sample
-   Approximation of Ridges and Umbilics on Triangulated Surface Meshes
    (new package)
    This package enables the approximation of differential features on
    triangulated surface meshes. Such curvature related features are
    lines: ridges or crests, and points: umbilics.
-   Envelopes of Curves in 2D (new package)
    This package contains two sets of functions that construct the lower
    and upper envelope diagram for a given range of bounded or unbounded
    curves.
-   Envelopes of Surfaces in 3D (new package)
    This package contains two sets of functions that construct the lower
    and upper envelope diagram for a given range of bounded or unbounded
    surfaces. The envelope diagram is realized as a 2D arrangement.
-   Minkowski Sums in 2D (new package)
    This package contains functions for computing planar Minkowski sums
    of two closed polygons, and for a polygon and a disc (an operation
    also known as offsetting or dilating a polygon). The package also
    contains an efficient approximation algorithm for the offset
    computation, which provides a guaranteed approximation bound while
    significantly expediting the running times w.r.t. the exact
    computation procedure.
-   Surface Mesh Parametrization
    Added Jacobi and SSOR preconditioners to OpenNL solver, which makes
    it much faster and more stable.
-   2D Arrangements
    -   Added support for unbounded curves.
    -   Added a traits class that supports bounded and unbounded linear
        objects, namely lines, rays and line segments.
    -   Added traits classes that handle circular arcs based on the
        circular kernel.
    -   Added a traits class that supports Bezier curves.
    -   Enhanced the traits class that supports rational functions to
        handle unbounded (as well as bounded) arcs
    -   Added a free function called decompose() that produces the
        symbolic vertical decomposition of a given arrangement,
        performing a batched vertical ray-shooting query from all
        arrangement vertices.
    -   Fixed a memory leak in the sweep-line code.
    -   Fixed a bug in computing the minor axis of non-degenerate
        hyperbolas.
-   Boolean Set Operations
    -   Added the DCEL as a default template parameter to the
        General\_polygon\_set\_2 and Polygon\_set\_2 classes. This
        allows users to extend the DCEL of the underlying arrangement.
    -   Added a function template called connect\_holes() that connects
        the holes in a given polygon with holes, turning it into a
        sequence of points, where the holes are connected to the outer
        boundary using zero-width passages.
    -   Added a non-const function member to General\_polygon\_set\_2
        that obtains the underlying arrangement.
-   2D and 3D Triangulations
    -   The constructors and insert member functions which take an
        iterator range perform spatial sorting in order to speed up the
        insertion.
-   Optimal Distances
    -   Polytope\_distance\_d: has support for homogeneous points;
        bugfix in fast exact version.
-   Bounding Volumes
    -   Min\_annulus\_d has support for homogeneous points; bugfix in
        fast exact version.

### Support Library

-   CGAL and the Boost Graph Library (BGL) (new package)
    This package provides the glue layer for several CGAL data
    structures such that they become models of the BGL graph concept.
-   Spatial Sorting (new package)
    This package allows to sort points and other objects along a Hilbert
    curve which can improve the performance of algorithms like
    triangulations. It is used by the constructors of the triangulation
    package which have an iterator range of points as argument.
-   Linear and Quadratic Programming Solver (new package)
    This package contains algorithms for minimizing linear and convex
    quadratic functions over polyhedral domains, described by linear
    equations and inequalities.

Release 3.2.1
-------------

Release date: July 2006

This is a bug fix release

### Number Types

-   Fix MP\_Float constructor which crashed for some values.

### Kernel

-   Rename Bool to avoid a clash with a macro in X11 headers.

### Arrangement

-   Derived the Arr\_segment\_traits\_2 Arrangement\_2 traits class from
    the parameterized Kernel. This allows the use of this traits class
    in an extended range of applications that require kernel objects and
    operations on these objects beyond the ones required by the
    Arrangement\_2 class itself.
-   Fixed a compilation bug in the code that handles overlay of
    arrangements instantiated with different DCEL classes.
-   Fixed a couple of bugs in the implementation of the Trapezoidal RIC
    point-location strategy

### Triangulation, Alpha Shapes

-   Qualify calls to filter\_iterator with "CGAL::" to avoid overload
    ambiguities with Boost's filter\_iterator.

### Surface Mesher

-   Fixed a bug in iterators of the class template
    Surface\_mesh\_complex\_2\_in\_triangulation\_3

### Surface Mesh Parametrisation

-   Updated the precompiled taucs lib

### Kinetic Data Structures

-   Fixed problems caused by old versions of gcc being confused by
    operator! and operator int()
-   Added point removal support to the Active\_objects\_vector

Release 3.2
-----------

Release date: May 2006

Version 3.2 differs from version 3.1 in the platforms that are supported
and in functionality. There have also been a number of bug fixes for
this release.

The following platforms are no longer supported:

-   SunPro CC versions 5.4 and 5.5 on Solaris
-   SGI Mips Pro

For Visual C++ the installation scripts choose the multi-threaded
dynamically linked runtime (/MD). Before it was the single-threaded
static runtime (/ML).

### Installation

-   The install tool tries to find third party libraries at "standard"
    locations.
-   Installers for Apple, Windows, and rpms.

### Manuals

-   User and Reference manual pages of a package are in the same chapter

### Kernels

-   2D Circular Kernel (new package)
    This package is an extension of the linear CGAL Kernel. It offers
    functionalities on circles, circular arcs and line segments in the
    plane.

### Data Structures and Algorithms

-   2D Regularized Boolean Set-Operations (new package)
    This package consists of the implementation of Boolean
    set-operations on point sets bounded by weakly x-monotone curves in
    2-dimensional Euclidean space. In particular, it contains the
    implementation of regularized Boolean set-operations, intersection
    predicates, and point containment predicates.
-   2D Straight Skeleton and Polygon Offsetting (new package)
    This package implements an algorithm to construct a halfedge data
    structure representing the straight skeleton in the interior of 2D
    polygons with holes and an algorithm to construct inward offset
    polygons at any offset distance given a straight skeleton.
-   2D Voronoi Diagram Adaptor (new package)
    This package provides an adaptor that adapts a 2-dimensional
    triangulated Delaunay graph to the corresponding Voronoi diagram,
    represented as a doubly connected edge list (DCEL) data structure.
    The adaptor has the ability to automatically eliminate, in a
    consistent manner, degenerate features of the Voronoi diagram, that
    are artifacts of the requirement that Delaunay graphs should be
    triangulated even in degenerate configurations. Depending on the
    type of operations that the underlying Delaunay graph supports, the
    adaptor allows for the incremental or dynamic construction of
    Voronoi diagrams and can support point location queries.
-   3D Surface Mesher (new package)
    This package provides functions to generate surface meshes that
    interpolate smooth surfaces. The meshing algorithm is based on
    Delaunay refinement and provides some guarantees on the resulting
    mesh: the user is able to control the size and shape of the mesh
    elements and the accuracy of the surface approximation. There is no
    restriction on the topology and number of components of input
    surfaces. The surface mesher may also be used for non smooth
    surfaces but without guarantee.

    Currently, implementations are provided for implicit surfaces
    described as the zero level set of some function and surfaces
    described as a gray level set in a three-dimensional image.

-   3D Surface Subdivision Methods (new package)
    Subdivision methods recursively refine a control mesh and generate
    points approximating the limit surface. This package consists of
    four popular subdivision methods and their refinement hosts.
    Supported subdivision methods include Catmull-Clark, Loop, Doo-Sabin
    and sqrt(3) subdivisions. Their respective refinement hosts are PQQ,
    PTQ, DQQ and sqrt(3) refinements. Variations of those methods can be
    easily extended by substituting the geometry computation of the
    refinement host.
-   Planar Parameterization of Triangulated Surface Meshes (new
    package)
    Parameterizing a surface amounts to finding a one-to-one mapping
    from a suitable domain to the surface. In this package, we focus on
    triangulated surfaces that are homeomorphic to a disk and on
    piecewise linear mappings into a planar domain. This package
    implements some of the state-of-the-art surface mesh
    parameterization methods, such as least squares conformal maps,
    discrete conformal map, discrete authalic parameterization, Floater
    mean value coordinates or Tutte barycentric mapping.
-   Principal Component Analysis (new package)
    This package provides functions to compute global information on
    the shape of a set of 2D or 3D objects such as points. It provides
    the computation of axis-aligned bounding boxes, centroids of point
    sets, barycenters of weighted point sets, as well as linear least
    squares fitting for point sets in 2D, and point sets as well as
    triangle sets in 3D.
-   2D Placement of Streamlines (new package)
    Visualizing vector fields is important for many application domains.
    A good way to do it is to generate streamlines that describe the
    flow behavior. This package implements the "Farthest Point Seeding"
    algorithm for placing streamlines in 2D vector fields. It generates
    a list of streamlines corresponding to an input flow using a
    specified separating distance. The algorithm uses a Delaunay
    triangulation to model objects and address different queries, and
    relies on choosing the centers of the biggest empty circles to start
    the integration of the streamlines.
-   Kinetic Data Structures (new package)
    Kinetic data structures allow combinatorial structures to be
    maintained as the primitives move. The package provides
    implementations of kinetic data structures for Delaunay
    triangulations in two and three dimensions, sorting of points in one
    dimension and regular triangulations in three dimensions. The
    package supports exact or inexact operations on primitives which
    move along polynomial trajectories.
-   Kinetic Framework (new package)
    Kinetic data structures allow combinatorial geometric structures to
    be maintained as the primitives move. The package provides a
    framework to ease implementing and debugging kinetic data
    structures. The package supports exact or inexact operations on
    primitives which move along polynomial trajectories.
-   Smallest Enclosing Ellipsoid (new package)
    This algorithm is new in the chapter Geometric Optimization.
-   2D Arrangement (major revision)
    This package can be used to construct, maintain, alter, and display
    arrangements in the plane. Once an arrangement is constructed, the
    package can be used to obtain results of various queries on the
    arrangement, such as point location. The package also includes
    generic implementations of two algorithmic frameworks, that are,
    computing the zone of an arrangement, and line-sweeping the plane,
    the arrangements is embedded on.

    Arrangements and arrangement components can also be extended to
    store additional data. An important extension stores the
    construction history of the arrangement, such that it is possible to
    obtain the originating curve of an arrangement subcurve.

-   Geometric Optimization (major revision)
    The underlying QP solver which is the foundation for several
    algorithms in the Geometric Optimization chapter has been completely
    rewritten.
-   3D Triangulation (new functionality)
    Regular\_triangulation\_3 now offers vertex removal.

Release 3.1
-----------

Release date: December 2004

Version 3.1 differs from version 3.0 in the platforms that are supported
and in functionality. There have also been a number of bug fixes for
this release.

Additional supported platforms:

-   MS Visual C++, version 7.3. and 8.0
-   Intel 8.0
-   SunPro CC versions 5.4 and 5.5 on Solaris
-   GNU g++ versions 3.4 on Linux, Solaris, Irix, cygwin, FreeBSD, and
    MacOS X
-   Darwin (MacOS X) and IA64/Linux support.

The following platforms are no longer supported:

-   MS Visual C++, version 7.0

The following functionality has been added or changed:


### All

-   The CORE 1.7 library for exact real arithmetic.
-   Updated GMP to 4.1.3.
-   Added Mpfr a library for multiple-precision floating-point
    computations with exact rounding.
-   Added Boost 1.32.0 (only include files).

### Installation

-   new option --disable-shared to omit building libCGAL.so.

### Manuals

-   Merged all major manuals in one multi-part manual, which provides
    now cross-links between the CGAL Kernel, the CGAL Basic Library, and
    the CGAL Support Library HTML manuals.
-   Improved layout.

### Kernels

-   Improved efficiency of filtered kernels.
-   More predicates and constructions.

### Basic Library

-   2D Segment Voronoi Diagram (new package)
    A data structure for Voronoi diagrams of segments in the plane under
    the Euclidean metric. The Voronoi edges are arcs of straight lines
    and parabolas. The algorithm provided in this package is
    incremental.
-   2D Conforming Triangulations and Meshes (new package)
    An implementation of Shewchuk's algorithm to construct conforming
    triangulations and 2D meshes.
-   3D Boolean Operations on Nef Polyhedra (new package)
    A new class (Nef\_polyhedron\_3) representing 3D Nef polyhedra, a
    boundary representation for cell-complexes bounded by halfspaces
    that supports boolean operations and topological operations in full
    generality including unbounded cells, mixed dimensional cells (e.g.,
    isolated vertices and antennas). Nef polyhedra distinguish between
    open and closed sets and can represent non-manifold geometry.
-   2D and Surface Function Interpolation (new package)
    This package implements different methods for scattered data
    interpolation: Given measures of a function on a set of discrete
    data points, the task is to interpolate this function on an
    arbitrary query point. The package further offers functions for
    natural neighbor interpolation.
-   Planar Nef polyhedra embedded on the sphere (new package)
    A new class (Nef\_polyhedron\_S2) designed and supported mainly to
    represent sphere neighborhoods around vertices of the three-
    dimensional Nef polyhedra.
-   Box\_intersection\_d (new package)
    A new efficient algorithm for finding all intersecting pairs for
    large numbers of iso-oriented boxes, i.e., typically these will be
    bounding boxes of more complicated geometries. Useful for (self-)
    intersection tests of surfaces etc.
-   2D Snap Rounding (new package)
    Snap Rounding is a well known method for converting
    arbitrary-precision arrangements of segments into a fixed-precision
    representation. In the study of robust geometric computing, it can
    be classified as a finite precision approximation technique.
    Iterated Snap Roundingis a modification of Snap Rounding in which
    each vertex is at least half-the-width-of-a-pixel away from any
    non-incident edge. This package supports both methods.
-   3D Triangulations
    -   Triangulation\_3: added operator==(),removed push\_back() and
        copy\_triangulation().
    -   Delaunay\_3 : added nearest\_vertex(), move\_point(),
        vertices\_in\_conflict().
    -   Regular\_3 : added filtered traits class, and
        nearest\_power\_vertex().
-   Planar\_map and Arrangement\_2
    -   The interface of the two traits functions that compute the
        intersection of two given curves changed. The functions
        nearest\_intersection\_to\_right() and
        nearest\_intersection\_to\_left() return an object of type
        CGAL::Object that represents either an empty intersection, a
        point, or an overlapping subcurve.
    -   Requirements to define two binary tags were added to the traits
        concept of the Planar\_map as follows: *Has\_left\_category* -
        indicates whether the functions
        curves\_compare\_y\_at\_x\_left() and
        nearest\_intersection\_to\_left() are implemented in the traits
        model. *Has\_reflect\_category* - indicates whether the
        functions point\_reflect\_in\_x\_and\_y() and
        curve\_reflect\_in\_x\_and\_y() are implemented in the traits
        model. They can be used as an alternative to the two function in
        the previous item.
    -   A new constructor of the Segment\_cached\_2 type that represents
        a segment in the Arr\_segment\_cached\_traits\_2 traits class
        was introduced. The new constructor accepts the segment
        endpoints as well as the coefficients of the underlying line.
    -   A new version of the conic-arc traits, based on CORE version 1.7
        was introduced. This new traits class makes use of CORE's
        rootOf() operator to compute the intersection points in the
        arrangement, making its code much simpler and more elegant than
        the previous version. In addition, new constructors for conic
        arcs are provided. The new traits class usually performs about
        30% faster than the version included in CGAL 3.0
    -   The traits class that handles continuous piecewise linear
        curves, namely Arr\_polyline\_traits\_2, was rewritten. The new
        class is parametrized with a traits class that handles segments,
        say Segment\_traits. The polyline curve defined within the
        Arr\_polyline\_traits\_2 class is implemented as a vector of
        segments of type Segment\_traits::Curve\_2.
    -   A meta traits class, namely Arr\_curve\_data\_traits\_2, that
        extends the curve type of the planar-map with arbitrary
        additional data was introduced. It should be instantiated with a
        regular traits-class and a class that contains all extraneous
        data associated with a curve.
    -   The class that represents the trapezoidal-decomposition point
        location strategy was renamed to
        Pm\_trapezoid\_ric\_point\_location.
    -   The Arrangement demo was rewritten. It covers many more
        features, has a much better graphical user interface, and comes
        with online documentation.
    -   Few bugs in the sweep-line module related to overlapping
        vertical segments were fixed. This module is used by the
        aggregate insert method that inserts a collection of curves at
        once.
-   Triangulation\_2
    -   added a filtered trait class in the regular triangulation
    -   added split and join operations in the triangulation data
        structure class
-   Alpha\_shapes\_3
    -   major changes in the implementation of the class
        Alpha\_shapes\_3.
    -   New implementation results in a true GENERAL mode allowing null
        and negative alpha-values. It also fixed the edges
        classification bug and introduces a classification of vertices.
-   Min\_ellipse\_2
    -   made access to approximate double representation public
    -   fixed bugs in conversion to double representation
    -   added `is_circle()` method
    -   minor performance improvements
-   Min\_sphere\_of\_spheres\_d:
    -   The models
        `Min_sphere_of_spheres_d_traits_2<K,FT,UseSqrt,Algorithm>`,
        `Min_sphere_of_spheres_d_traits_3<K,FT,UseSqrt,Algorithm>`, and
        `Min_sphere_of_spheres_d_traits_d<K,FT,Dim,UseSqrt,Algorithm>`
        of concept `MinSphereOfSpheresTraits` now represent a sphere as
        a `std::pair<Point,Radius>` (and not any more as a
        `CGAL::Weighted_point<Point,Weight>`)
    -   Internal code cleanup; in particular, implementation details
        don't pollute the namespace CGAL anymore
-   Polyhedron\_3
    -   New Tutorial on CGAL Polyhedron for Subdivision Algorithms with
        interactive demo viewer in source code available.
    -   Added example program for efficient self-intersection test. -
        Added small helper functions, such as vertex\_degree,
        facet\_degree, edge\_flip, and is\_closed.
-   Apollonius Graph (Voronoi of Circles)
    -   Reduced memory requirements by approximately a factor of two.

Release 3.0.1
-------------

Release date: February 2004

This is a bug-fix release. No new features have been added in 3.0.1.
Here is the list of bug-fixes.

### Polyhedral Surface

-   Fixed wrong include files for output support. Added example.

### Planar\_map

-   Fixed the so called "Walk-along-a-line" point-location strategy to
    correctly handle a degenerate case.

### 2D Triangulation

-   added missing figure in html doc
-   in Line\_face\_circulator\_2.h:
    Fixed changes made to support handles with a typedef to iterator.
    The fix concerns operator== and !=.

### Alpha\_shapes\_3

-   fixed classify member function for edges.

### Number types

-   Lazy\_exact\_nt:
    -   added the possibility to select the relative precision of
        `to_double()` (by default 1e-5). This should fix reports that
        some circumcenters computations have poor coordinates, e.g.
        nan).
    -   when exact computation is triggered, the interval is recomputed,
        this should speed up some kinds of computations.
-   `to_interval(Quotient<MP_Float>)`: avoid spurious overflows.

### Kernel

-   missing acknowledgment in the manual and minor clarification of
    `intersection()` documentation.

Release 3.0
-----------

Release date: October 2003

Version 3.0 differs from version 2.4 in the platforms that are supported
and in functionality. There have also been a number of bug fixes for
this release.

The license has been changed to either the LGPL (GNU Lesser General
Public License v2.1) or the QPL (Q Public License v1.0) depending on
each package. So CGAL remains free of use for you, if your usage meets
the criteria of these licenses, otherwise, a commercial license has to
be purchased from GeometryFactory.

Additional supported platforms:

-   MS Visual C++, version 7.1.
-   SunPro CC versions 5.4 and 5.5 on Solaris
-   GNU g++ versions 3.2 and 3.3 on Linux, Solaris, Irix, cygwin, and
    FreeBSD.
-   MipsPRO CC 7.30 and 7.40 with both the n32 and n64 ABIs.

The following platforms are no longer supported:

-   MS Visual C++, version 6.
-   GNU g++ 2.95.2 (2.95.3 is still supported)
-   Kai C++ and Borland C++, all versions

The following functionality has been added or changed:

**All**

-   The CORE library for exact computations is now distributed as part
    of CGAL as well.

### Kernels

-   3 typedefs have been added to ease the choice of a robust and fast
    kernel:
    -   Exact\_predicates\_inexact\_constructions\_kernel
    -   Exact\_predicates\_exact\_constructions\_kernel
    -   Exact\_predicates\_exact\_constructions\_kernel\_with\_sqrt
-   Progress has been made towards the complete adaptability and
    extensibility of our kernels.
-   New faster Triangle\_3 intersection test routines.
    *(see Erratum)*
-   Added a Kernel concept archetype to check that generic algorithms
    don't use more functionality than they should.
-   A few more miscellaneous functions.

### Basic Library

-   2D Apollonius Graph (new package)
    Algorithms for computing the Apollonius graph in two dimensions. The
    Apollonius graph is the dual of the Apollonius diagram, also known
    as the additively weighted Voronoi diagram. The latter can be
    thought of as the Voronoi diagram of a set of circles under the
    Euclidean metric, and it is a generalization of the standard Voronoi
    diagram for points. The algorithms provided are dynamic.
-   dD Min Sphere of Spheres (new package)
    Algorithms to compute the smallest enclosing sphere of a given set
    of spheres in R<sup>d</sup>. The package provides an algorithm with
    maximal expected running time *O(2<sup>O(d)</sup> n)* and a fast and
    robust heuristic (for dimension less than 30).
-   Spatial Searching (new package)
    Provides exact and approximate distance browsing in a set of points
    in *d*-dimensional space using implementations of algorithms
    supporting:
    -   both nearest and furthest neighbor searching
    -   both exact and approximate searching
    -   (approximate) range searching
    -   (approximate) *k*-nearest and *k*-furthest neighbor searching
    -   (approximate) incremental nearest and incremental furthest
        neighbor searching
    -   query items representing points and spatial objects.
-   **Kd-tree**
    this package is deprecated, its documentation is removed. It is
    replaced by the Spatial Searching package.
-   Largest\_empty\_rectangle\_2
    Given a set of points P in the plane, the class
    Largest\_empty\_iso\_rectangle\_2 is a data structure that maintains
    an iso-rectangle with the largest area among all iso-rectangles that
    are inside a given iso-rectangle bounding box, and that do not
    contain any point of the point set P.
-   2D Triangulation and 3D Triangulation
    -   The classes Triangulation\_data\_structure\_2 (and 3), which
        implements the data structure for 2D triangulation class, now
        makes use of CGAL::Compact\_container (see Support Library
        section below).
    -   The triangulation classes use a Rebind mechanism to provide the
        full flexibility on Vertex and Face base classes. This means
        that it is possible for the user to derive its own Face of
        Vertex base class, adding a functionality that makes use of
        types defined by the triangulation data structure like
        Face\_handle or Vertex\_handle.
    -   New classes Triangulation\_vertex\_base\_with\_info\_2 (and 3)
        and Triangulation\_face\_base\_with\_info\_2 (and 3) to make
        easier the customization of base classes in most cases.
-   2D Triangulation
    -   Regular triangulation provides an easy access to hidden points.
    -   The Triangulation\_hierarchy\_2, which provide an efficient
        location data structure, can now be used with any 2D
        triangulation class plugged in (including Regular
        triangulations).
-   3D Triangulation
    -   faster vertex removal function in Delaunay\_triangulation\_3.
    -   Delaunay\_triangulation\_3 is now independent of the order of
        insertions of the points (in case of degenerate cosphericity).
    -   Regular\_triangulation\_3 now hides vertices (and updates
        itself) when inserting a coinciding point with greater weight.
        This required a new predicate.
    -   deprecated functions: copy\_triangulation(), push\_back(),
        set\_number\_of\_vertices().
    -   Triangulation\_3 now gives non-const access to the data
        structure.
-   Interval Skip List (new package)
    An interval skip list is a data structure for finding all intervals
    that contain a point, and for stabbing queries, that is for
    answering the question whether a given point is contained in an
    interval or not.
-   Planar Maps and Arrangements
    The changes concern mainly the traits classes.
    1.  New traits hierarchy and interface: The set of requirements was
        made sound and complete. A couple of requirements were
        eliminated, few others were redefined, and some were renamed. A
        hierarchy of three traits classes for the Planar\_map\_2,
        Planar\_map\_with\_intersections\_2, and Arrangement\_2 types
        was established to include only the necessary requirements at
        each level. It was determined that for the aggregate insertion-
        operation based on a sweep-line algorithm only a subset of the
        requirements is needed. Preconditions were added where
        appropriate to tighten the requirements further.

        The following functions have been renamed:

        -   point\_is\_same() renamed to point\_equal()
        -   curve\_is\_same() renamed to curve\_equal()
        -   curve\_is\_in\_x\_range() renamed to point\_in\_x\_range()
        -   curve\_compare\_at\_x() renamed to
            curves\_compare\_y\_at\_x() Furthermore, a precondition has
            been added that the reference point is in the x-range of
            both curves.
        -   curve\_compare\_at\_x\_right() renamed to
            curves\_compare\_y\_at\_x\_to\_right(). Furthermore, a
            precondition has been added that both curves are equal at
            the reference point and defined to its right.
        -   curve\_compare\_at\_x\_left() renamed to
            curves\_compare\_y\_at\_x\_to\_left(). Furthermore, a
            precondition has been added that both curves are equal at
            the reference point and defined to its right.
        -   curve\_get\_point\_status() renamed to
            curve\_compare\_y\_at\_x(). Furthermore, a precondition has
            been added that the point is in the x-range of the curve.
            Consequently, the function now returns a Comparison\_result
            (instead of a special enum).
        -   make\_x\_monotone() renamed to curve\_make\_x\_monotone()
            See more details below.
        -   curve\_flip() renamed to curve\_opposite()

        The following functions have been removed:

        -   curve\_is\_between\_cw()
        -   point\_to\_left()
        -   point\_to\_right()
        -   is\_x\_monotone()
        -   point\_reflect\_in\_x\_and\_y()
        -   curve\_reflect\_in\_x\_and\_y()
        -   do\_intersect\_to\_right()
        -   do\_intersect\_to\_left()

        Most functions, are required by the PlanarMapTraits\_2 concept,
        except for the make\_x\_monotone(),
        nearest\_intersection\_to\_right(),
        nearest\_intersection\_to\_left(), curves\_overlap() and
        curve\_opposite(). PlanarMapWithIntersectionsTraits\_2 requires
        all these functions, except curve\_opposite(), needed only by
        the ArrangementTraits\_2 concept.

        Furthermore, the two functions curve\_compare\_at\_x\_left() and
        nearest\_intersection\_to\_left() can be omitted, if the two
        functions point\_reflect\_in\_x() and curve\_reflect\_in\_x()
        are implemented. Reflection can be avoided, if the two \_left
        functions are supplied.

    2.  The type X\_curve\_2 of the PlanarMapWithIntersectionsTraits\_2
        concept was renamed to X\_monotone\_curve\_2, and the
        distinction between this type and the Curve\_2 type was made
        firm. The method is\_x\_monotone() of the
        PlanarMapWithIntersectionsTraits\_2 concept was removed. The
        related method curve\_make\_x\_monotone() is now called for each
        input curve of type Curve\_2 when curves are inserted into a
        Planar\_map\_with\_intersections\_2 to subdivide the input curve
        into x-monotone sub-curves (and in case the curve is already
        x-monotone, this function is responsible for casting it to an
        x-monotone curve).
    3.  New and improved traits classes:
    4.  Conic traits - Arr\_conic\_traits\_2 Support finite segments of
        ellipses, hyperbolas and parabolas, as well as line segments.
        The traits require an exact real number- type, such as
        leda\_real or CORE::Expr.
    5.  Segment cached traits - Arr\_segment\_cached\_traits\_2 This
        class uses an improved representation for segments that helps
        avoiding cascaded computations, thus achieving faster running
        times. To work properly, an exact rational number-type should be
        used.
    6.  Polyline traits - Arr\_polyline\_traits\_2 The polyline traits
        class has been reimplemented to work in a more efficient,
        generic manner. The new class replaces the obsolete
        Arr\_polyline\_traits class. It is parameterized with a segment
        traits class.
    7.  Hyperbola and segment traits - Arr\_hyper\_segment\_traits\_2
        Supports line segments and segments of canonical hyperbolas.
        This is the type of curves that arise when projecting segments
        in three-space rotationally around a line onto a plane
        containing the line. Such projections are often useful in
        CAD/CAM problems.
    8.  Removed old traits class:
        -   The models of the PlanarMapWithIntersectionsTraits\_2
            concept below became obsolete, as the new conic traits,
            namely Arr\_conic\_traits\_2, supports the same
            functionality and is much more efficient.
            -   Arr\_circles\_real\_traits
            -   Arr\_segment\_circle\_traits
        -   The segment traits class and the new polyline traits class
            were reimplemented using standard CGAL-kernel calls. This
            essentially eliminated the corresponding leda traits
            classes, namely:
            -   Pm\_leda\_segment\_traits\_2
            -   Arr\_leda\_segment\_traits\_2
            -   Arr\_leda\_polyline\_traits

            With the use of the Leda\_rat\_kernel new external package
            the same functionality can be achieved with less overhead
            and more efficiency.
    9.  Sweep Line
        -   The Sweep\_line\_2 package was reimplemented. As a
            consequence it is much more efficient, its traits is tighter
            (namely neither the two \_left nor the reflection functions
            are required), and its interface has changed a bit.
            1.  The following global functions have been removed:
                -   sweep\_to\_produce\_subcurves\_2()
                -   sweep\_to\_produce\_points\_2()
                -   sweep\_to\_construct\_planar\_map\_2()

                Instead, the public methods of the Sweep\_line\_2 class
                listed below were introduced:
                -   get\_subcurves() - Given a container of curves, this
                    function returns a list of curves that are created
                    by intersecting the input curves.
                -   get\_intersection\_points() - Given a range of
                    curves, this function returns a list of points that
                    are the intersection points of the curves.
                -   get\_intersecting\_curves() - Given a range of
                    curves, this function returns an iterator to the
                    beginning of a range that contains the list of
                    curves for each intersection point between any two
                    curves in the specified range.

            2.  It is possible to construct a planar map with
                intersections (or an arrangement) by inserting a range
                of curves into an empty map. This will invoke the
                sweep-line process to construct the map more
                efficiently.
        -   New interface functions to the
            Planar\_map\_with\_intersections\_2 class. The
            Planar\_map\_with\_intersections\_2 class maintains a planar
            map of input curves that possibly intersect each other and
            are not necessarily x-monotone. If an input curve, or a set
            of input curves, are known to be x-monotone and pairwise
            disjoint, the new functions below can be used to insert them
            into the map efficiently.

-   Polyhedral Surface
    -   The old design that was deprecated since CGAL 2.3 has been
        removed.
    -   Class `Polyhedron_incremental_builder_3`:
        -   Renamed local enum `ABSOLUTE` to `ABSOLUTE_INDEXING`, and
            `RELATIVE` to `RELATIVE_INDEXING` to avoid conflicts with
            similarly named macros of another library.
        -   Changed member functions `add_vertex()`, `begin_facet()`,
            and `end_facet()` to return useful handles.
        -   Added `test_facet()` to check facets for validity before
            adding them.
        -   Added `vertex( size_t i)` to return `Vertex_handle` for
            index `i`.
-   Halfedge Data Structure
    -   The old design that was deprecated since CGAL 2.3 has been
        removed.

### Support Library

-   New container class Compact\_container, which (roughly) provides the
    flexibility of std::list, with the memory compactness of
    std::vector.
-   Geomview\_stream: added a function gv.draw\_triangles(InputIterator
    begin, InputIterator end) which draws a set of triangles much more
    quickly than one by one.
-   Number types:
    -   number types are now required to provide a function:
        std::pair&lt;double, double&gt; to\_interval(const NT &).
    -   number types are now required to provide mixed operators with
        "int".
    -   CLN support removed.
    -   faster square() for MP\_Float.
    -   added Gmp\_q.
-   Qt\_widget:
    -   New classes:
        -   Qt\_help\_window: provides a simple way to show some helpful
            information about a demo as an HTML page.
        -   Qt\_widget\_history: provides basic functionality to
            manipulate intervals of Qt\_widget class. The current
            visible area of Qt\_widget is mapped to an interval. Each
            interval could be stored in the Qt\_widget\_history object.
            So you can use this object to navigate in history. It is
            mostly used by Qt\_widget\_standard\_toolbar.
    -   Changes:
        -   Qt\_widget\_standard\_toolbar: is derived from QToolBar
            class, so pay attention to modify your code, if you used
            this class. Some public methods were introduced to control
            the history object that the toolbar use to navigate.
        -   the icons are now part of libCGALQt.
    -   Deprecated members of Qt\_widget:
        -   add\_to\_history(), clear\_history(), back(), forth(): use
            forward(), back() and clear\_history() of the
            Qt\_widget\_standard\_toolbar instead.
        -   custom\_redraw(): use redraw\_on\_back() and
            redraw\_on\_front() instead.
    -   Optimizations: the output operators of the following classes
        have been optimized:
        -   CGAL::Segment\_2 (now tests for intersection with the
            drawing area)
        -   CGAL::Triangle\_2 (now tests for intersection with the
            drawing area)
        -   CGAL::Triangulation\_2 (is optimized for faster display on
            zooming)

### Erratum in the Kernel manual

-   Intersection test routines

    The documentation of CGAL::do\_intersect should mention, for the 3D
    case:
    Also, in three-dimensional space *Type1* can be

    -   either *Plane\_3&lt;Kernel&gt;*
    -   or *Triangle\_3&lt;Kernel&gt;*

    and *Type2* any of

    -   *Plane\_3&lt;Kernel&gt;*
    -   *Line\_3&lt;Kernel&gt;*
    -   *Ray\_3&lt;Kernel&gt;*
    -   *Segment\_3&lt;Kernel&gt;*
    -   *Triangle\_3&lt;Kernel&gt;*

    In the same way, for *Kernel::DoIntersect\_3*:
    for all pairs *Type1* and *Type2*, where the type *Type1* is

    -   either *Kernel::Plane\_3*
    -   or *Kernel::Triangle\_3*

    and *Type2* can be any of the following:

    -   *Kernel::Plane\_3*
    -   *Kernel::Line\_3*
    -   *Kernel::Ray\_3*
    -   *Kernel::Segment\_3*
    -   *Kernel::Triangle\_3*

    Philippe Guigue (INRIA Sophia-Antipolis) should be mentioned as one
    of the authors.

Release 2.4
-----------

Release date: May 2002

Version 2.4 differs from version 2.3 in the platforms that are supported
and in functionality. There have also been a number of bug fixes for
this release.

Additional supported platforms:

-   Microsoft Visual C++, version 7.
-   SunPro 5.3 (with patch 111685-05) on Solaris
-   g++ 3.1 on Linux and Solaris

The following functionality has been added or changed:


### Kernels

-   Point\_d has been removed from the 2D and 3D kernels. This type is
    now available from the d-dimensional kernel only.

### Basic Library

-   2D Polygon Partitioning
    Traits requirements for optimal partitioning have been changed
    slightly.

-   2D Sweep line
    A new package that implements a sweep-line algorithm to compute
    arrangements of curves for different families of curves, which are
    not necessarily line segments (e.g., it also works for circular
    arcs). The resulting output can be the list of vertex points, the
    resulting subcurves or a planar map.

-   Planar Maps and Arrangements
    -   New quicker insertion functions of Planar\_map\_2 for cases
        where more precomputed information is available regarding the
        position of the inserted curve in the map.
    -   New query function for planar maps that determines whether a
        given point is within a given face of the planar map.
    -   New iterator over edges of planar maps in addition to the
        existing iterator over halfedges.
    -   New copy constructor and assignment operator for arrangements.



-   Polyhedral Surface
    -   new design introduced with release 2.3 now supported by VC7
        compiler
    -   Extended functionality of Polyhedron\_incremental\_builder:
        absolute indexing allows one to add new surfaces to existing
        ones.



-   2D Triangulation
    -   There is a new triangulation data structure replacing the two
        previous ones. This new data structure is coherent with the 3d
        triangulation data structure and offer the advantages of both
        previous ones. Backward compatibility is ensured and this change
        is transparent for the user of triangulation classes.
    -   Constrained and Delaunay constrained triangulations are now able
        to handle intersecting input constraints. The behavior of
        constrained triangulations with respect to intersection of input
        constraints can be customized using an intersection tag.
    -   A new class Constrained\_triangulation\_plus offers a
        constrained hierarchy on top of a constrained triangulations.
        This additional data structure describes the subdivision of the
        original constraints into edges of the triangulations.



-   3D Triangulation
    -   Running time improved by a better and more compact management of
        memory allocation
    -   Various improvements and small functionalities added:
        -   Triangulation\_3&lt;GT,Tds&gt;::triangle() returns a
            triangle oriented towards the outside of the cell c for
            facet (c,i)
        -   New function insert(Point, Locate\_type, Cell\_handle, int,
            int) which avoids the location step.
        -   New function to get access to cells in conflict in a
            Delaunay insertion : find\_conflicts() and
            insert\_in\_hole()
        -   New function TDS::delete\_cells(begin, end).
        -   New functions : degree(v), reorient(),
            remove\_decrease\_dimension(), remove\_from\_simplex().
    -   Changes of interface:
        -   vertices and cells are the same for the triangulation data
            structure and the geometric triangulation
        -   the triangulation data structure uses Vertex\_handle (resp
            Cell\_handle) instead of Vertex\* (resp Cell\*).
        -   incident\_cells() and incident\_vertices() are templated by
            output iterators
        -   changes in the iterators and circulators interface:
            -   Iterators and circulators are convertible to handles
                automatically, no need to call "-&gt;handle()" anymore.
            -   Vertex\_iterator split into All\_vertices\_iterator and
                Finite\_vertices\_iterator (and similar for cells...).
            -   TDS::Edge/Facet iterators now support operator-&gt;.



-   2D Search structures
    Additional range search operations taking a predicate functor have
    been added

### Support Library

-   Qt\_widget
    -   We have added a new class for visualization of 2D CGAL objects.
        It is derived from Trolltech's Qt class QWidget and privdes a
        used to scale and pan.
    -   Some demos were developed for the following packages: 2D Alpha
        shapes, 2D Convex Hull, Largest empty 2D rectangle, Maximum
        k-gon, Minimum ellipse, Minimum 2D quadrilateral, 2D polygon
        partitioning 2D regular and constrained triangulation.
    -   Tutorials are available to help users get used to Qt\_widget



-   Timer
    Fixed Timer class (for user process time) to have no wrap-around
    anymore on Posix-compliant systems.

The following functionality is no longer supported:

-   Planar maps of infinite curves (the so-called planar map
    bounding-box).

Bugs in the following packages have been fixed: 3D Convex hull, 2D
Polygon partition, simple polygon generator

Also attempts have been made to assure compatibility with the upcoming
LEDA release that introduces the leda namespace.

### Known problems

-   2D Nef Polyhedra contains a memory leak. Memory problems are also
    the likely cause of occasional run-time errors on some platforms.
-   The d-dimensional convex hull computation produces run-time errors
    on some platforms because of memory management bugs.
-   The new Halfedge Data Structure design introduced with release 2.3
    does not work on VC6. See the release notes in the manual for more
    information.
-   The following deficiencies relate to planar maps, planar maps of
    intersecting curves (pmwx), arrangements and sweep line.
    -   On KCC, Borland and SunPro we guarantee neither compilation nor
        correct execution for all of the packages above.
    -   On VC6 and VC7 we guarantee neither compilation nor correct
        execution of the sweep line package.
    -   On CC (on Irix 6.5) the trapezoidal decomposition point location
        strategy is problematic when used with planar maps, pmwx, or
        arrangements (mind that this is the default for planar maps).
    -   On CC (on Irix 6.5) sweep line with polyline traits does not
        compile (mind that the so-called leda polyline traits does
        compile).
    -   On g++ (on Irix 6.5) the segment-circle
        (Arr\_segment\_circle\_traits\_2) traits does not compile for
        either of the above packages.

Release 2.3
-----------

Release date: August 2001

Version 2.3 differs from version 2.2 in the platforms that are supported
and in functionality.

Additional supported platform:

-   Gnu g++ 3.0 on Solaris and Linux

The following functionality has been added:


### Kernels

-   The 2D and 3D kernels now serve as models of the new kernel concept
    described in the recent paper, "An Adaptable and Extensible Geometry
    Kernel" by Susan Hert, Micheal Hoffmann, Lutz Kettner, Sylvain Pion,
    and Michael Seel to be presented at WAE 2001 (and soon available as
    a technical report). This new kernel is completely compatible with
    the previous design but is more flexible in that it allows geometric
    predicates as well as objects to be easily exchanged and adapted
    individually to users' needs.
-   A new kernel called `Simple_homogeneous` is available. It is
    equivalent to `Homogeneous` but without reference-counted objects.
-   A new kernel called `Filtered_kernel` is available that allows one
    to build kernel traits classes that use exact and efficient
    predicates.
-   There are two classes, `Cartesian_converter` and
    `Homogeneous_converter` that allows one to convert objects between
    different Cartesian and homogeneous kernels, respectively.
-   A new d-dimensional kernel, `Kernel_d` is available. It provides
    diverse kernel objects, predicates and constructions in d dimensions
    with two representations based on the kernel families `Cartesean_d`
    and `Homogeneous_d`

### Basic Library

Almost all packages in the basic library have been adapted to the new
kernel design to realize the flexibility this design makes possible. In
several packages, this means that the traits class requirements have
changed to conform to the function objects offered in the kernels so the
kernels themselves can be used as traits classes in many instances.
-   2D Convex Hull
    The traits requirements have changed slightly to bring them in line
    with the CGAL kernels.
-   3D Convex Hull
    -   The function `convex_hull_3` now uses a new implementation of
        the quickhull algorithm and no longer requires LEDA.
    -   A new `convex_hull_incremental_3` function based on the new
        d-dimensional convex hull class is available for comparison
        purposes.


-   `Convex_hull_d, Delaunay_d`
    Two new application classes offering the calculation of
    d-dimensional convex hulls and delaunay triangulations

-   Polygons and Polygon Operations
    -   The traits class requirements have been changed.
    -   The simplicity test has a completely new implementation.
    -   Properties like convexity, simplicity and area can now be cached
        by polygons. You need to set a flag to select this behavior.




-   Planar Nef Polyhedra
    A new class (`Nef_polyhedron_2`) representing planar Nef polyhedra =
    rectilinearly bounded points sets that are the result of binary and
    topological operations starting from halfplanes.

-   A new package offering functions to partition planar polygons into
    convex and y-monotone pieces is available.

-   Planar Maps and Arrangements
    -   A new class `Planar_map_with_intersections_2<Planar_map>` for
        planar maps of possibly intersecting, possibly non-x-monotone,
        possibly overlapping curves (like `Arrangement_2` but without
        the hierarchy tree).
    -   I/O utilities for planar maps and arrangements for textual and
        graphical streams. (It is possible to save and later reload
        built planar maps or arrangements.)
    -   New arrangement traits class for line segments and circular arcs
        (`Arr_segment_circle_traits<NT>`).
    -   New faster traits for polylines specialized for using the LEDA
        rational kernel (`Arr_leda_polylines_traits`). The LEDA traits
        for segments was also made faster.
    -   A new point location strategy
        (`Pm_simple_point_location<Planar_map>`).



-   Halfedge Data Structure

    The halfedge data structure has been completely revised. The new
    design is more in line with the STL naming scheme and it provides a
    safe and coherent type system throughout the whole design (no void\*
    pointers anymore), which allows for better extendibility. A user can
    add new incidences in the mesh easily. The new design also uses
    standard allocators with a new template parameter that has a
    suitable default.

    The old design is still available, but its use is deprecated, see
    the manual of deprecated packages for its documentation. Reported
    bugs in copying the halfedge data structure (and therefore also
    polyhedral surfaces) have been fixed in both designs. Copying a
    list-based representation is now based on hash maps instead of
    std::map and is therefore considerably faster.

-   Polyhedral Surface

    The polyhedral surface has been rewritten to work with the new
    halfedge data structure design. The user level interface of the
    `CGAL::Polyhedron_3` class is almost backwards compatible with the
    previous class. The exceptions are the template parameter list,
    everything that relies on the flexibility of the underlying halfedge
    data structure, such as a self-written facet class, and that the
    distinction between supported normals and supported planes has been
    removed. Only planes are supported. See the manuals for suggestions
    how to handle normals instead of planes.

    More example programs are provided with polyhedral surfaces, for
    example, one about Euler operator and one computing a subdivision
    surface given a control mesh as input.

    The old design is still available for backwards compatibility and to
    support older compiler, such as MSVC++6.0. For the polyhedral
    surface, old and new design cannot be used simultaneously (they have
    identical include file names and class names). The include files
    select automatically the old design for MSVC++6.0 and the new design
    otherwise. This automatism can be overwritten by defining
    appropriate macros before the include files. The old design is
    selected with the `CGAL_USE_POLYHEDRON_DESIGN_ONE` macro. The new
    design is selected with the `CGAL_USE_POLYHEDRON_DESIGN_TWO`
    macro.

-   2D Triangulation
    -   The geometric traits class requirements have been changed to
        conform to the new CGAL kernels. CGAL kernel classes can be used
        as traits classes for all 2D triangulations except for regular
        triangulations.
    -   Additional functionality:
        -   dual method for regular triangulations (to build a power
            diagram)
        -   unified names and signatures for various "find\_conflicts()"
            member functions in Delaunay and constrained Delaunay
            triangulation.
        -   As an alternative to the simple insert() member function,
            insertion of points in those triangulation can be performed
            using the combination of find\_conflicts() and star\_hole()
            which eventually allows the user to keep track of deleted
            faces.
    -   More demos and examples


-   3D Triangulation
    -   Major improvements
        -   A new class `Triangulation_hierarchy_3` that allows a faster
            point location, and thus construction of the Delaunay
            triangulation
        -   A new method for removing a vertex from a Delaunay
            triangulation that solves all degenerate cases
        -   Running time of the usual location and insertion methods
            improved
    -   A bit more functionality, such as
        -   New geomview output
        -   dual methods in Delaunay triangulations to draw the Voronoi
            diagram
    -   More demos and examples
    -   Changes in interface
        -   Traits classes requirements have been modified
        -   The kernel can be used directly as a traits class (except
            for regular triangulation)
        -   insert methods in `Triangulation_data_structure` have a new
            interface


-   A new class (`Alpha_shapes_3`) that computes Alpha shapes of point
    sets in 3D is available.

-   The traits requirements for matrix search and minimum quadrilaterals
    have been changed to bring them in line with the CGAL kernels.

-   Point\_set\_2
    -   now independent of LEDA; based on the CGAL Delaunay
        triangulation
    -   traits class requirements adapted to new kernel concept.
    -   function template versions of the provided query operations are
        available

### Support Library

-   Number types:
    -   `Lazy_exact_nt<NT>` is a new number type wrapper to speed up
        exact number types.
    -   `MP_Float` is a new multiprecision floating point number type.
        It can do exact additions, subtractions and multiplications over
        floating point values.
-   `In_place_list` has a new third template parameter (with a suitable
    default) for an STL-compliant allocator.
-   `Unique_hash_map` is a new support class.
-   `Union_find` is a new support class.
-   `Geomview_stream` :
    -   Geomview version 1.8.1 is now required.
    -   no need to have a `~/.geomview` file anymore.
    -   new output operators for triangulations.
    -   new output operators for `Ray_2`, `Line_2`, `Ray_3`, `Line_3`,
        `Sphere_3`.
    -   various new manipulators.
-   Window stream In cooperation with Algorithmic Solutions, GmBH
    (distributors of the LEDA library), we can now offer a visualization
    package downloadable in binary form that supports visualization on a
    ported version of the LEDA window lib.

Release 2.2
-----------

Release date: October 2000

Version 2.2 differs from version 2.1 in the platforms that are supported
and in functionality.

Additional supported platforms:

-   the KAI compiler (4.0) on Solaris 5.8
-   Borland C++ (5.5)

The following functionality has been added:

-   There is a new, non-reference-counted kernel, Simple\_cartesian.
    Because reference counting is not used, and thus coordinates are
    stored within a class, debugging is easier using this kernel. This
    kernel can also be faster in some cases than the reference-counted
    Cartesian kernel.
-   New optimization algorithms
    -   Min\_annulus\_d - Algorithm for computing the smallest enclosing
        annulus of points in arbitrary dimension
    -   Polytope\_distance\_d - Algorithm for computing the (squared)
        distance between two convex polytopes in arbitrary dimension
    -   Width\_3 - Algorithm for computing the (squared) width of points
        sets in three dimensions
-   2D Triangulations
    -   There are now two triangulation data structures available in
        CGAL. The new one uses a list to store the faces and allows one
        to represent two-dimensional triangulations embedded in three
        spaces as well as planar triangulations.
    -   The triangulation hierarchy which allows fast location query is
        now available.
-   Infinite objects can now be included in planar maps.
-   Removal as well as insertions of vertices for 3D Delaunay
    triangulations is now possible.
-   A generator for \`\`random'' simple polygons is now available.
-   In directory demo/Robustness, programs that demonstrate typical
    robustness problems in geometric computing are presented along with
    the solutions to these problems that CGAL provides.

The following functionality has been removed:

-   The binary operations on polygons (union, intersection ...) have
    been removed. Those operations were not documented in the previous
    release (2.1). Arrangements can often be used as a substitute.

Release 2.1
-----------

Release date: January 2000

Version 2.1 differs from version 2.0 in the platforms that are supported
and in functionality.

Supported platforms:

-   the newest gnu compiler (2.95.2) on Sun, SGI, Linux and Windows.
-   the Microsoft Visual C++ compiler, version 6.
-   the mips CC compiler version 7.3 under Irix.

Support for the old g++ compiler (2.8) and for mips CC 7.2 has been
dropped.

The following functionality has been added:
-   Alpha shapes and weighted alpha shapes in 2D. Alpha shapes are a
    generalization of the convex hull of a point set.
-   Arrangements in 2D. Arrangements are related to and based on planar
    maps. The major difference between the two is that curves are
    allowed to intersect in the case of arrangements.
-   Extensions to triangulations in 2D. Constrained triangulations are
    now dynamic: they support insertions of new constraint as well as
    removal of existing constraints. There are also constrained Delaunay
    triangulations.
-   Triangulations in 3D were added, both Delaunay triangulations and
    regular triangulations.
-   Min\_quadrilateral optimizations have been added. These are
    algorithms to compute the minimum enclosing rectangle/parallelogram
    (arbitrary orientation) and the minimum enclosing strip of a convex
    point set.
-   2d Point\_set is a package for 2d range search operations, Delaunay
    triangulation, nearest neighbor queries. This package works only if
    LEDA is installed.
-   Support for GeoWin visualization library. This also depends on LEDA.
-   Support for using the CLN number type together with CGAL.

Release 2.0
-----------

Release date: June 1999

The main difference from release 1.2 is the introduction of namespaces
-- namespace `std` for code from the standard library and namespace
`CGAL` for the CGAL library.

Release 1.2
-----------

Release date: January 1999

Additions to release 1.1 include:

-   topological map
-   planar map overlay
-   regular and constrained triangulations

Release 1.1
-----------

Release date: July 1998

Additions to release 1.0 include:

-   3D intersections
-   kD points
-   3D convex hull
-   kD smallest enclosing sphere

Release 1.0
-----------

Release date: April 1998

Additions to release 0.9 include:

-   Polyhedral surfaces
-   Halfedge Data Structure
-   Planar maps

Release 0.9
-----------

Release date: June 1997

Initial (beta) release of the CGAL library.
