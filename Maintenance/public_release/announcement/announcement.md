%The CGAL Open Source Project is pleased to announce the release 5.3 of CGAL, the Computational Geometry Algorithms Library.

Besides fixes and general enhancement to existing packages, the following has changed since CGAL 5.2:

### [General changes](https://doc.cgal.org/5.3/Manual/general_intro.html)

-   The support for the compiled version of CGAL is dropped. Only the header-only version is supported.

-   On Windows, the type used for `Exact_rational`, in `Epick` and indirectly (through `Lazy_exact_nt`)
   `Epeck` may now be `boost::multiprecision::mpq_rational`, as has been the case on other platforms
   for several releases. This depends on various options and is added to a list that includes
   `mpq_class`, `CGAL::Gmpq`, `leda_rational` and `CGAL::Quotient<CGAL::MP_Float>`.

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
    which can be used to determine whehter a closed path on a surface mesh can be continously
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
