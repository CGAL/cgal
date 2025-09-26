The CGAL Open Source Project is pleased to announce the release 6.1 Beta 1 of CGAL, the Computational Geometry Algorithms Library.

CGAL version 6.1 Beta 1 is a public testing release. It should provide a solid ground to report bugs that need to be tackled before the release of the final version of CGAL 6.1 in Sept 2025.

Besides fixes and general enhancement to existing packages, the following has changed since CGAL 6.0:

### General Changes

- The new list of supported compilers is:
  - Visual C++ 15.9, 16.10, 17.0 (from Visual Studio 2017, 2019 and 2022) or later
  - Gnu g++ 12.2.0 or later (on Linux)
  - LLVM Clang version 20.1.6 or later (on Linux)
  - Apple Clang compiler versions 12.0.5 and 12.0.5 (on macOS)
- The minimal supported version of Boost is now 1.74.0.

### [3D Constrained Triangulations](https://doc.cgal.org/6.1/Manual/packages.html#PkgConstrainedTriangulation3) (new package)

- This package implements the construction of a 3D Constrained Delaunay triangulation. This triangulation is a generalization of a 3D Delaunay Triangulation which conforms to the set of faces of a 3D piecewise linear complex (PLC), ensuring that these faces are part of the triangulation. As not all PLCs are tetrahedralizable, the algorithm may insert Steiner points to construct the constrained triangulation. The main entry point is the function [`CGAL::make_conforming_constrained_Delaunay_triangulation_3()`](https://doc.cgal.org/6.1/Constrained_triangulation_3/group__PkgConstrainedTriangulation3FunctionsPolygonSoupOrMesh.html).

### [3D Isosurfacing](https://doc.cgal.org/6.1/Manual/packages.html#PkgIsosurfacing3) (new package)

- This package provides algorithms to extract isosurfaces from scalar fields. The algorithms provided in this first version include Marching Cubes, Topologically Correct Marching Cubes, and Dual Contouring. The algorithm is generic with respect to the scalar field representation (implicit function, discrete values, …) and the discretization data structure (Cartesian grid, octree, …). The output is an indexed face set that stores an isosurface in the form of a surface mesh.

### [dD Fréchet Distance](https://doc.cgal.org/6.1/Manual/packages.html#PkgFrechetDistance) (new package)

- This package provides functions for computing the Fréchet distance of polylines in any dimension under the Euclidean metric.

### [2D Triangulations on Hyperbolic Surfaces](https://doc.cgal.org/6.1/Manual/packages.html#PkgHyperbolicSurfaceTriangulation2) (new package)

- This package enables building and handling triangulations of closed orientable hyperbolic surfaces. It offers functions for the generation of the triangulation from a convex fundamental domain, the Delaunay flip algorithm, and the construction of a portion of the lift of the triangulation in the Poincaré disk. A method is offered that generates such domains in genus two.

  See also the associated [news entry](https://www.cgal.org/2025/06/24/triangulations-on-hyperbolic-surfaces/).

### [Polygon Repair](https://doc.cgal.org/6.1/Manual/packages.html#PkgPolygonRepair)

- Added the [non-zero rule](https://doc.cgal.org/6.1/Polygon_repair/structCGAL_1_1Polygon__repair_1_1Non__zero__rule.html) (areas with non-zero winding number are kept), as well as two functions to compute the conservative inner and outer hull of similar polygons:
  - [`CGAL::Polygon_repair::join()`](https://doc.cgal.org/6.1/Polygon_repair/group__PkgPolygonRepairFunctions.html#gad5b959666d952392c0e3b8d4b1b1847a)
  - [`CGAL::Polygon_repair::intersect()`](https://doc.cgal.org/6.1/Polygon_repair/group__PkgPolygonRepairFunctions.html#ga780e31115643e3d0b406349b56c9f3d5)

  See also the associated [news entry](https://www.cgal.org/2025/05/22/Polygon_repair/).

### [Polygon Mesh Processing](https://doc.cgal.org/6.1/Manual/packages.html#PkgPolygonMeshProcessing)

- Added the parameter `apply_iterative_snap_rounding` to the function [`CGAL::Polygon_mesh_processing::autorefine_triangle_soup()`](https://doc.cgal.org/6.1/Polygon_mesh_processing/group__PMP__corefinement__grp.html#gaf7747d676c459d9e5da9b13be7d12bb5). When set to `true`, the coordinates are rounded to fit in double and may perform additional subdivisions to ensure the output is free of self-intersections. See also the associated [news entry](https://www.cgal.org/2025/06/13/autorefine-and-snap/).
- Added the function [`CGAL::Polygon_mesh_processing::approximated_centroidal_Voronoi_diagram_remeshing()`](https://doc.cgal.org/6.1/Polygon_mesh_processing/group__PkgPolygonMeshProcessingRef.html#gaed23e63b12c7fe8268927d17b4d379f1) to remesh triangle meshes. This remeshing algorithm uses clustering on polygonal meshes as to approximate a Centroidal Voronoi Diagram construction, and can move vertices as to recover sharp features and corners. See also the associated [news entry](https://www.cgal.org/2025/05/22/Surface_remeshing/).
- New implementation of [`CGAL::Polygon_mesh_processing::clip()`](https://doc.cgal.org/6.1/Polygon_mesh_processing/group__PMP__corefinement__grp.html#ga2c73d3460872e601f84a03f58dd069ae) and [`CGAL::Polygon_mesh_processing::split()`](https://doc.cgal.org/6.1/Polygon_mesh_processing/group__PMP__corefinement__grp.html#ga6738052411a4d548a5b375f11b913924) with a plane as clipper that is much faster and is now able to handle non-triangulated surface meshes. See also the associated [news entry](https://www.cgal.org/2025/06/06/new_clip/).
- Added the function [`CGAL::Polygon_mesh_processing::refine_with_plane()`](https://doc.cgal.org/6.1/Polygon_mesh_processing/group__PMP__corefinement__grp.html#gacb9d68fa4dea8fd03ec53b56a16d6fc6), which enables users to refine a mesh with its intersection with a plane.
- Added a function in the [visitor of the corefinement based methods](https://doc.cgal.org/6.1/Polygon_mesh_processing/classPMPCorefinementVisitor.html) to trace faces in the output meshes which correspond to coplanar faces of the input.
- Added the function [`CGAL::Polygon_mesh_processing::discrete_mean_curvature()`](https://doc.cgal.org/6.1/Polygon_mesh_processing/group__PMP__measure__grp.html#ga1a31fa9412b4643dc7202a54246db78b) and [`CGAL::Polygon_mesh_processing::discrete_Gaussian_curvature()`](https://doc.cgal.org/6.1/Polygon_mesh_processing/group__PMP__measure__grp.html#ga11a2d646d4636605d185653bff5bbbbb) to evaluate the discrete curvature at a vertex of a mesh.
- Added the function [`CGAL::Polygon_mesh_processing::angle_sum()`](https://doc.cgal.org/6.1/Polygon_mesh_processing/group__PMP__measure__grp.html#ga25d3836c21931610cf76b6716a06254c) to compute the sum of the angles around a vertex.

### [Point Set Processing](https://doc.cgal.org/6.1/Manual/packages.html#PkgPointSetProcessing3)

- Added [`CGAL::poisson_eliminate()`](https://doc.cgal.org/6.1/Point_set_processing_3/group__PkgPointSetProcessing3Algorithms.html#ga2d73d46ca766656a219bf5e6045b837a), which can be used to downsample a point cloud to a target size while providing Poisson disk property, i.e., a larger minimal distance between points.

### [CGAL and the Boost Graph Library (BGL)](https://doc.cgal.org/6.1/Manual/packages.html#PkgBGL)

- Added the function [`dijkstra_shortest_path()`](https://doc.cgal.org/6.1/BGL/group__PkgBGLTraversal.html#gaa4058482db0089886b84a8c6a341e528), which can be used to compute the geometrically shortest sequence of halfedges between two vertices.
- Added the function [`CGAL::Euler::remove_degree_2_vertex()`](https://doc.cgal.org/6.1/BGL/group__PkgBGLEulerOperations.html#gab3455663b7db4624529e54ae3ce2387c), which enables users to remove vertices which have exactly two incident edges.

### [2D Arrangements](https://doc.cgal.org/6.1/Manual/packages.html#PkgArrangementOnSurface2)

- **Breaking change**: Renamed the concept `AosApproximateTraits_2` to [`AosApproximatePointTraits_2`](https://doc.cgal.org/6.1/Arrangement_on_surface_2/classAosApproximatePointTraits__2.html) to make room for the new concept [`AosApproximateTraits_2`](https://doc.cgal.org/6.1/Arrangement_on_surface_2/classAosApproximateTraits__2.html). This concept requires the provision of a functor called `Approximate_2` that has an operator that approximates the coordinates of a point.
- **Breaking change**: The concept [`AosApproximateTraits_2`](https://doc.cgal.org/6.1/Arrangement_on_surface_2/classAosApproximateTraits__2.html) now refines the concept [`AosApproximatePointTraits_2`](https://doc.cgal.org/6.1/Arrangement_on_surface_2/classAosApproximatePointTraits__2.html) and requires the provision of a functor called `Approximate_2`. In addition to an operator that approximates the coordinates of a point, it also requires the provision of (i) an operator that approximates a points, and (ii) an operator that approximates a curve.
- Renamed the prefix of the names of all concepts in the Arrangement_on_surface_2 package from “Arrangement” to “Aos”.
- Introduced two traits decorators, namely [`CGAL::Arr_tracing_traits_2`](https://doc.cgal.org/6.1/Arrangement_on_surface_2/classCGAL_1_1Arr__tracing__traits__2.html) and [`CGAL::Arr_counting_traits_2`](https://doc.cgal.org/6.1/Arrangement_on_surface_2/classCGAL_1_1Arr__counting__traits__2.html), which can be used to extract debugging and informative metadata about the traits in use while a program is being executed.
- Fixed the Landmark point-location strategy so that it can be applied to arrangements on a sphere.
- Fixed a bug in the extensions of vertex and halfedge types of the DCEL when used to instantiate Arrangement_with_history_2 or similar arrangement classes that derive from Arrangement_2.
- Fixed `do_intersect()` of a 2D Arrangement with a curve.

### Triangulations

- All triangulations now offer the functions `point(Vertex_handle)` and `point(Simplex, int)`, which enables users to access the geometric position of a vertex and of the i-th vertex of a simplex of a triangulation.

### [Poisson Surface Reconstruction](https://doc.cgal.org/6.1/Manual/packages.html#PkgPoissonSurfaceReconstruction3)

- Added a new mesh domain `Poisson_mesh_domain_3` that integrates some optimizations from the deprecated 3D Surface Mesh Generation package.

### [2D Triangulations](https://doc.cgal.org/6.1/Manual/packages.html#PkgTriangulation2)

- **Breaking change**: In the class template [`Constrained_triangulation_plus_2`](https://doc.cgal.org/6.1/Triangulation_2/classCGAL_1_1Constrained__triangulation__plus__2.html), the value type of the range returned by [`subconstraints()`](https://doc.cgal.org/6.1/Triangulation_2/classCGAL_1_1Constrained__triangulation__plus__2.html#af25114a7e1675194367f8f9de9de90d2) has changed from `const std::pair<const Subconstraint, std::list<Context>*>` to `Subconstraint`. The old range type is now returned by a new function named `subconstraints_and_contexts()`.

### [3D Mesh Generation](https://doc.cgal.org/6.1/Manual/packages.html#PkgMesh3)

- Added two new meshing parameters that enable custom mesh initialization:
- [`initial_points_generator`](https://doc.cgal.org/6.1/Mesh_3/group__PkgMesh3Parameters.html#gaf53777b83f1b2f3e7d49275dbab6e46b): enables the user to specify a functor that generates initial points,
- [`initial_points`](https://doc.cgal.org/6.1/Mesh_3/group__PkgMesh3Parameters.html#gaf26d164d1845dcd66dc4861b6920b5ec): enables the user to specify a `Range` of initial points.
- Added a new meshing parameter [`surface_only`](https://doc.cgal.org/6.1/Mesh_3/group__PkgMesh3Parameters.html#gaa2618c09b6117d7caab12dccca16ee58), which can be used to improve performance when only surface mesh generation is sought.
- Added a new mesh domain [`Poisson_mesh_domain_3`](https://doc.cgal.org/6.1/Poisson_surface_reconstruction_3/classCGAL_1_1Poisson__mesh__domain__3.html), which should be used when generating a mesh from a Poisson surface obtained with the package [Poisson Surface Reconstruction](https://doc.cgal.org/6.1/Manual/packages.html#PkgPoissonSurfaceReconstruction3). This mesh domain re-integrates some optimizations for Poisson surface mesh generation that were lost when the package [3D Mesh Generation](https://doc.cgal.org/6.1/Manual/packages.html#PkgMesh3) had to be replaced instead of the deprecated package [3D Surface Mesh Generation](https://doc.cgal.org/latest/Manual/packages.html#PkgSurfaceMesher3).

### [3D Subdivision Methods](https://doc.cgal.org/6.1/Manual/packages.html#PkgSurfaceSubdivisionMethod3)

- Added a new named parameter for [`CGAL::Subdivision_method_3::Loop_subdivision()`](https://doc.cgal.org/6.1/Subdivision_method_3/group__PkgSurfaceSubdivisionMethod3Functions.html#gafa1e441c4e07eb06e1f6efecef7ff268) and [`CGAL::Subdivision_method_3::CatmullClark_subdivision()`](https://doc.cgal.org/6.1/Subdivision_method_3/group__PkgSurfaceSubdivisionMethod3Functions.html#ga8e6c8dd3c26d7a27c070b3a091684679), which enables users to subdivide a mesh without modifying its geometry.

### [Algebraic Kernel](https://doc.cgal.org/6.1/Manual/packages.html#PkgAlgebraicKernelD)

- **Breaking change**: Classes based on the RS Library are no longer provided.
