The CGAL Open Source Project is pleased to announce the release 5.4 of CGAL, the Computational Geometry Algorithms Library.

Besides fixes and general enhancement to existing packages, the following has changed since CGAL 5.3:

### [General changes](https://doc.cgal.org/5.4/Manual/general_intro.html)

-   Added the cmake target `CGAL::CGAL_Basic_viewer` to ease the compilation of programs using the basic viewer-based function `CGAL::draw()`. This target will define the macro and link with `CGAL_Qt5` target when linked with it.

-   The kernel providing exact constructions and exact predicates ([`CGAL::Exact_predicates_exact_constructions_kernel`](https://doc.cgal.org/5.4/Kernel_23/classCGAL_1_1Exact__predicates__exact__constructions__kernel.html)) is now thread-safe. See changes in `2D and 3D Linear Geometry Kernel` for more details.

-   The class `Geomview_stream` and all the dependent features have been removed from CGAL. Those features were actually no longer supported since CGAL-5.3 but it was not properly announced.

### [Shape Regularization](https://doc.cgal.org/5.4/Manual/packages.html#PkgShapeRegularization) (new package)

-   This package enables to regularize a set of segments and open or closed contours in 2D and a set of planes in 3D such that all input objects are rotated and aligned with respect to the user-specified conditions. In addition, it provides a global regularization framework that can be adjusted for the user needs and any type of geometric objects.

### [Weights](https://doc.cgal.org/5.4/Manual/packages.html#PkgWeights) (new package)

-   This package provides a simple and unified interface to different types of weights. In particular, it groups all weights into three category: analytic weights including all basic weights which can be computed analytically for a query point with respect to its local neighbors in 2D and 3D; barycentric weights, including all weights which can be computed for a query point with respect to the vertices of a planar polygon; and weighting regions, including all weights which are used to balance other weights.

### [2D Generalized Barycentric Coordinates](https://doc.cgal.org/5.4/Manual/packages.html#PkgBarycentricCoordinates2) (major changes)

-   **Breaking change**: The headers `Segment_coordinates_2.h` and `Triangle_coordinates_2.h` are renamed to `segment_coordinates_2.h` and `triangle_coordinates_2.h`.
-   The classes [`Segment_coordinates_2`](https://doc.cgal.org/5.4/Barycentric_coordinates_2/classCGAL_1_1Barycentric__coordinates_1_1Segment__coordinates__2.html) and [`Triangle_coordinates_2`](https://doc.cgal.org/5.4/Barycentric_coordinates_2/classCGAL_1_1Barycentric__coordinates_1_1Triangle__coordinates__2.html) are deprecated. The free functions [`compute_segment_coordinates_2()`](https://doc.cgal.org/5.4/Barycentric_coordinates_2/classCGAL_1_1Barycentric__coordinates_1_1Segment__coordinates__2.html#a134d363dccaeecb5621fa608fac76eaf) and [`compute_triangle_coordinates_2()`](https://doc.cgal.org/5.4/Barycentric_coordinates_2/classCGAL_1_1Barycentric__coordinates_1_1Triangle__coordinates__2.html#a958fee3ad9613d7bfa9d7a976aa3548f) are deprecated as well. Instead, the free functions [`segment_coordinates_2()`](https://doc.cgal.org/5.4/Barycentric_coordinates_2/group__PkgBarycentricCoordinates2RefFunctions.html#gab856ca68d37f58e6cdf74c8aac6f4245) and [`triangle_coordinates_2()`](https://doc.cgal.org/5.4/Barycentric_coordinates_2/group__PkgBarycentricCoordinates2RefFunctions.html#gaa378786f8996dbcefe7923ebb711e4dd) should be used.
-   The enums [`Query_point_location`](https://doc.cgal.org/5.4/Barycentric_coordinates_2/namespaceCGAL_1_1Barycentric__coordinates.html#aedeeb072a2024053a016afd15e591331) and [`Type_of_algorithm`](https://doc.cgal.org/5.4/Barycentric_coordinates_2/namespaceCGAL_1_1Barycentric__coordinates.html#a5e5682512438422f23d6080edc49c05b) are deprecated. Instead, the enum [`Computation_policy_2`](https://doc.cgal.org/5.4/Barycentric_coordinates_2/namespaceCGAL_1_1Barycentric__coordinates.html#a478bbcec416216b2274ee4b4e97b0e6c) should be used.
-   The classes [`Wachspress_2`](https://doc.cgal.org/5.4/Barycentric_coordinates_2/classCGAL_1_1Barycentric__coordinates_1_1Wachspress__2.html), [`Discrete_harmonic_2`](https://doc.cgal.org/5.4/Barycentric_coordinates_2/classCGAL_1_1Barycentric__coordinates_1_1Discrete__harmonic__2.html), [`Mean_value_2`](https://doc.cgal.org/5.4/Barycentric_coordinates_2/classCGAL_1_1Barycentric__coordinates_1_1Mean__value__2.html), and [`Generalized_barycentric_coordinates_2`](https://doc.cgal.org/5.4/Barycentric_coordinates_2/classCGAL_1_1Barycentric__coordinates_1_1Generalized__barycentric__coordinates__2.html) are deprecated. As consequence, the concept [`BarycentricCoordinates_2`](https://doc.cgal.org/5.4/Barycentric_coordinates_2/classCGAL_1_1Barycentric__coordinates_1_1BarycentricCoordinates__2.html) is deprecated as well. Instead, the classes [`Wachspress_coordinates_2`](https://doc.cgal.org/5.4/Barycentric_coordinates_2/classCGAL_1_1Barycentric__coordinates_1_1Wachspress__coordinates__2.html), [`Discrete_harmonic_coordinates_2`](https://doc.cgal.org/5.4/Barycentric_coordinates_2/classCGAL_1_1Barycentric__coordinates_1_1Discrete__harmonic__coordinates__2.html), and [`Mean_value_coordinates_2`](https://doc.cgal.org/5.4/Barycentric_coordinates_2/classCGAL_1_1Barycentric__coordinates_1_1Mean__value__coordinates__2.html) should be used.
-   Added the class [`Harmonic_coordinates_2`](https://doc.cgal.org/5.4/Barycentric_coordinates_2/classCGAL_1_1Barycentric__coordinates_1_1Harmonic__coordinates__2.html) to compute approximate harmonic coordinates in 2D. These coordinates satisfy all properties of barycentric coordinates inside any simple polygon.
-   Added a new concept [`DiscretizedDomain_2`](https://doc.cgal.org/5.4/Barycentric_coordinates_2/classCGAL_1_1Barycentric__coordinates_1_1DiscretizedDomain__2.html) and a model of this concept called [`Delaunay_domain_2`](https://doc.cgal.org/5.4/Barycentric_coordinates_2/classCGAL_1_1Barycentric__coordinates_1_1Delaunay__domain__2.html), which is based on the [Mesh 2](https://doc.cgal.org/5.4/Manual/packages.html#PkgMesh2) package. A model of this concept is required to use [`Harmonic_coordinates_2`](https://doc.cgal.org/5.4/Barycentric_coordinates_2/classCGAL_1_1Barycentric__coordinates_1_1Harmonic__coordinates__2.html).
-   Added free functions to compute Wachspress, discrete harmonic, and mean value coordinates.
-   All free functions and classes are now using ranges and property maps.

### [2D and 3D Linear Geometry Kernel](https://doc.cgal.org/5.4/Manual/packages.html#PkgKernel23)

-   Most operations on [`CGAL::Exact_predicates_exact_constructions_kernel`](https://doc.cgal.org/5.4/Kernel_23/classCGAL_1_1Exact__predicates__exact__constructions__kernel.html) objects are now thread-safe if [`CGAL::Exact_rational`](https://doc.cgal.org/5.4/Number_types/group__nt__cgal.html#ga0849ff44771b19582218ebdfa5614f64) is [`mpq_class`](https://doc.cgal.org/5.3/Number_types/classmpq__class.html) (from `GMPXX`), `boost::multiprecision::mpq_rational` or [`CGAL::Quotient<CGAL::MP_Float>`](https://doc.cgal.org/5.3/Number_types/classCGAL_1_1MP__Float.html). The objects are not atomic though, so the usual restrictions on avoiding race conditions apply. For users who do not use threads, this can be disabled with `CGAL_HAS_NO_THREADS`.

-   Added documentation for the class [`Projection_traits_3`](https://doc.cgal.org/5.4/Kernel_23/classCGAL_1_1Projection__traits__3.html), which enables the use of 2D algorithms on the projections of 3D data onto an arbitrary plane.

-   Added `construct_centroid_2_object()` and `compute_determinant_2_object()` in [`Projection_traits_xy_3`](https://doc.cgal.org/5.4/Kernel_23/classCGAL_1_1Projection__traits__xy__3.html), [`Projection_traits_xz_3`](https://doc.cgal.org/5.4/Kernel_23/classCGAL_1_1Projection__traits__xz__3.html), and [`Projection_traits_yz_3`](https://doc.cgal.org/5.4/Kernel_23/classCGAL_1_1Projection__traits__yz__3.html) classes.

-   Added the functor [`NonZeroCoordinateIndex_3`](https://doc.cgal.org/5.4/Kernel_23/classKernel_1_1NonZeroCoordinateIndex__3.html) to the concept [`Kernel`](https://doc.cgal.org/5.4/Kernel_23/classKernel.html) with `int operator()(Vector_3)` which returns the index of any coordinate of the vector different from zero, or `-1`.

### [dD Kernel](https://doc.cgal.org/5.4/Manual/packages.html#PkgKernelD)

-   Most operations on [`Epeck_d`](https://doc.cgal.org/5.4/Kernel_d/structCGAL_1_1Epeck__d.html) objects are now thread-safe, see 2D and 3D Linear Geometry Kernel for details.

### [2D Arrangements](https://doc.cgal.org/5.4/Manual/packages.html#PkgArrangementOnSurface2)

-   **Breaking Change:** The traits function objects `Compare_x_at_limit_2` and `Compare_x_near_limit_2` are renamed to `Compare_x_on_boundary_2` and `Compare_x_near_boundary_2`, respectively.

-   A [new hierarchy of traits concepts](https://doc.cgal.org/5.4/Arrangement_on_surface_2/group__PkgArrangementOnSurface2Concepts.html) has been introduced. It captures all the valid combinations of boundary conditions for the 4 boundary sides of the parameter space. The 4 boundaries are Bottom, Top, Left, and Right. Each boundary side can be either contracted, identified, close, open, or oblivious. Not all possible combinations are valid. If one side is identified then the other must be as well. Two adjacent sides cannot be contracted.

-   A new geometric traits, [`Arr_geodesic_arc_on_sphere_traits_2`](https://doc.cgal.org/5.4/Arrangement_on_surface_2/classCGAL_1_1Arr__geodesic__arc__on__sphere__traits__2.html) has been introduced. It handles arcs of great circles embedded on the unit sphere.

### [2D Regularized Boolean Set-Operations](https://doc.cgal.org/5.4/Manual/packages.html#PkgBooleanSetOperations2)

-   Added an extra parameter (`UsePolylines`) to all free functions ( [`complement()`](https://doc.cgal.org/5.4/Boolean_set_operations_2/group__boolean__complement.html), [`do_intersect()`](https://doc.cgal.org/5.4/Boolean_set_operations_2/group__boolean__do__intersect.html), [`intersection()`](https://doc.cgal.org/5.4/Boolean_set_operations_2/group__boolean__intersection.html), [`join()`](https://doc.cgal.org/5.4/Boolean_set_operations_2/group__boolean__join.html), [`difference()`](https://doc.cgal.org/5.4/Boolean_set_operations_2/group__boolean__difference.html), [`symmetric_difference()`](https://doc.cgal.org/5.4/Boolean_set_operations_2/group__boolean__symmetric__difference.html), and [`oriented_side`](https://doc.cgal.org/5.4/Boolean_set_operations_2/group__boolean__oriented__side.html)) to control whether to use `Arr_polyline_traits_2` as default traits. It is the new default as it provides better performances in general.

### [3D Mesh Generation](https://doc.cgal.org/5.4/Manual/packages.html#PkgMesh3)

-   Added support of weighted images for an improved quality of meshes generated from labeled images, along with a function [`CGAL::Mesh_3::generate_label_weights()`](https://doc.cgal.org/5.4/Mesh_3/namespaceCGAL_1_1Mesh__3.html#ae5914bf77180ff8948c08046154ee727) to generate the weights.

### [Polygon Mesh Processing](https://doc.cgal.org/5.4/Manual/packages.html#PkgPolygonMeshProcessing)

-   Added the function [`CGAL::Polygon_mesh_processing::match_faces()`](https://doc.cgal.org/5.4/Polygon_mesh_processing/group__measure__grp.html#ga10f7cd81645bafe936ac5eb4e58e67ef), which, given two polygon meshes, identifies their common faces as well as faces present in only either of them.

-   Added the functions: [`CGAL::Polygon_mesh_processing::bounded_error_Hausdorff_distance()`](https://doc.cgal.org/5.4/Polygon_mesh_processing/group__PMP__distance__grp.html#ga6d4ecea831c33ac10eec42b5021fc183) that computes an estimate of the one-sided Hausdorff distance between two triangle meshes which is bounded by a user-specified error bound; [`CGAL::Polygon_mesh_processing::bounded_error_symmetric_Hausdorff_distance()`](https://doc.cgal.org/5.4/Polygon_mesh_processing/group__PMP__distance__grp.html#ga9a7a682b5d9523135c8502e72117dffd) that computes an estimate of the symmetric Hausdorff distance bounded by a user-specified error bound; and [`CGAL::Polygon_mesh_processing::is_Hausdorff_distance_larger()`](https://doc.cgal.org/5.4/Polygon_mesh_processing/group__PMP__distance__grp.html#gab19e751107025a443e86baa9763aebf3) that returns `true` if the bounded-error Hausdorff distance between two meshes is larger than the user-specified max distance.

-   Added the functions [`CGAL::Polygon_mesh_processing::squared_edge_length()`](https://doc.cgal.org/5.4/Polygon_mesh_processing/group__measure__grp.html#ga30fa03722cd7aa599f6dcb115f54fec5) and [`CGAL::Polygon_mesh_processing::squared_face_area()`](https://doc.cgal.org/5.4/Polygon_mesh_processing/group__measure__grp.html#ga6eda3738815fd678df225f79ccfc3e03), which, compared to [`CGAL::Polygon_mesh_processing::edge_length()`](https://doc.cgal.org/5.4/Polygon_mesh_processing/group__measure__grp.html#gae1674775d9fecada7f25710f425cff3a) and [`CGAL::Polygon_mesh_processing::face_area()`](https://doc.cgal.org/5.4/Polygon_mesh_processing/group__measure__grp.html#ga6a1d7a825c09490b1e6613295343482b), enable avoiding square-root operations.

-   Added more functions in the [visitor of the corefinement based methods](https://doc.cgal.org/5.4/Polygon_mesh_processing/classPMPCorefinementVisitor.html) to track all vertex creations.

-   Added an option to [`CGAL::Polygon_mesh_processing::self_intersections()`](https://doc.cgal.org/5.4/Polygon_mesh_processing/group__PMP__intersection__grp.html#gaf19c80ec12cbff7ebe9e69453f1d40b8) to report only a limited number of intersections (`maximum_number()`).

### [The Heat Method](https://doc.cgal.org/5.4/Manual/packages.html#PkgHeatMethod)

-   **Breaking change**: Added the functor `Compute_squared_length_3` providing `operator(const Vector_3& v)`, which computes the squared length of `v`, to the [`HeatMethodTraits_3`](https://doc.cgal.org/5.4/Heat_method_3/classHeatMethodTraits__3.html) concept.

### [Point Set Processing](https://doc.cgal.org/5.4/Manual/packages.html#PkgPointSetProcessing3)

-   Added support for [`libpointmatcher::GenericDescriptorOutlierFilter`](https://github.com/ethz-asl/libpointmatcher) that enables providing a map from a point to a weight associated with this point.

### [Shape Detection](https://doc.cgal.org/5.4/Manual/packages.html#PkgShapeDetection)

-   Added new shapes to the Region Growing algorithm on a point set: circles in 2D, spheres in 3D, and cylinders in 3D.

### [CGAL and Solvers](https://doc.cgal.org/5.4/Manual/packages.html#PkgSolverInterface)

-   Added support for the [OSQP solver](https://osqp.org/). This solver enables to efficiently compute the convex Quadratic Programming (QP) problems arising in the context of several packages.
