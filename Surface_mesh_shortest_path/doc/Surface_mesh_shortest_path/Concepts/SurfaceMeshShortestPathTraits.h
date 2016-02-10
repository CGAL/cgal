/*!
\ingroup PkgSurfaceMeshShortestPathConcepts

\cgalConcept

The concept `SurfaceMeshShortestPathTraits` describes the types,
predicates, and constructions required by the traits class parameter of
`CGAL::Surface_mesh_shortest_path`.

\cgalHasModel `CGAL::Surface_mesh_shortest_path_traits<K,P>`

*/

class SurfaceMeshShortestPathTraits
{
public:

/// \name Types
/// @{

  /// A model of the concept `FaceListGraph`
  typedef unspecified_type Triangle_mesh;

  /// A model of the concept `FieldWithSqrt` or a model of both `Field` and `RealEmbeddable`.
  typedef unspecified_type FT;

  /// The 2-dimensional point type
  typedef unspecified_type Point_2;

  /// The 2-dimensional vector type
  typedef unspecified_type Vector_2;

  /// The 2-dimensional ray type
  typedef unspecified_type Ray_2;

  /// The 2-dimensional line type
  typedef unspecified_type Line_2;

  /// The 2-dimensional segment type
  typedef unspecified_type Segment_2;

  /// The 2-dimensional triangle type
  typedef unspecified_type Triangle_2;

  /// An ordered triple to specify barycentric coordinates in triangles.
  typedef unspecified_type Barycentric_coordinate;

  /// The 3-dimensional point type
  typedef unspecified_type Point_3;

  /// The 3-dimensional vector type
  typedef unspecified_type Vector_3;

  /// The 3-dimensional ray type
  typedef unspecified_type Ray_3;

  /// The 3-dimensional triangle type
  typedef unspecified_type Triangle_3;

/// @}

/// \name Constructions
/// @{

  /*!
  Function object type that provides
  `Point_2 operator()(FT x, FT y)` and
  `Point_2 operator()(CGAL::ORIGIN)` to construct points with
  cartesian coordinates `(x,y)` and `(0,0)` respectively.
  */
  typedef unspecified_type Construct_point_2;

  /*!
  Function object type that provides
  `Vector_2 operator()(Point_2 a, Point_2 b)`
  to construct the vector `b - a`
  */
  typedef unspecified_type Construct_vector_2;

  /*!
  Function object type that provides
  `Ray_2 operator()(Point_2 a, Point_2 b)`
  to construct the ray originating from `a` and passing through `b`
  */
  typedef unspecified_type Construct_ray_2;

  /*!
  Function object type that provides
  `Line_2 operator()(Segment_2 s)` and
  `Line_2 operator()(Ray_2 r)`
  to construct the supporting line to `s` and `r` respectively.
  */
  typedef unspecified_type Construct_line_2;

  /*!
  Function object type that provides
  `Segment_2 operator()(Point_2 a, Point_2 b)`
  to construct the segment `(a,b)`
  */
  typedef unspecified_type Construct_segment_2;

  /*!
  Function object type that provides
  `Triangle_2 operator()(Point_2 a, Point_2 b, Point_2 c)`
  to construct the triangle `(a,b,c)`
  */
  typedef unspecified_type Construct_triangle_2;

  /*!
  Function object type that provides
  `Point_2 operator()(Segment_2 s, int i)`
  to return the source or target of `s`, depending on whether `i` is 0 or 1 respectively, and
  `Point_2 operator()(Triangle_2 t, int i)`
  to return the `i`th vertex of `t`.
  */
  typedef unspecified_type Construct_vertex_2;

  /*!
  Function object type that provides
  `Point_2 operator()(Segment_2 s)`
  to return the source of `s`.
  */
  typedef unspecified_type Construct_source_2;

  /*!
  Function object type that provides
  `Point_2 operator()(Segment_2 s)`
  to return the target of `s`.
  */
  typedef unspecified_type Construct_target_2;

  /*!
  Function object type that provides
  `Point_2 operator()(Point_2 p1, FT w1, Point_2 p2, FT w2)`
  to compute the barycenter of the points `p1` and `p2` with corresponding weights `w1` and `w2`, and
  `Point_2 operator()(Point_2 p1, FT w1, Point_2 p2, FT w2, Point_2 p3, FT w3)`
  to compute the barycenter of the points `p1`, `p2`, and `p3` with corresponding weights `w1`, `w2`, and `w3`.
  */
  typedef unspecified_type Construct_barycenter_2;

  /*!
  Function object type that provides
  `FT operator()(A obj1, B obj2)`
  to compute the squared distance between `obj1` and `obj2`, `A` and `B` can be any one of
  `Point_2`,
  `Segment_2`,
  type objects.
  */
  typedef unspecified_type Compute_squared_distance_2;

  /*!
  Function object type.
  Must provide
  `CGAL::cpp11::result_of<Intersect_2(A,B)>::%type operator()(A obj1, B obj2)`
  to compute the intersection between `obj1` and `obj2`, where `A` and `B` can be any type amongst
  `Line_2`, `Ray_2`, `Segment_2`.
  */
  typedef unspecified_type Intersect_2;

  /*!
  Function object type that provides
  `Point_2 operator()(Ray_2 r, int i)`
  to construct a point on `r`, such that `i = 0` is the source, and any `i > 0` is not the source.
  */
  typedef unspecified_type Construct_point_on_2;

  /*!
  Function object type that provides
  `Vector_3 operator()(Point_3 a, Point_3 b)`
  to construct the vector `b - a`
  */
  typedef unspecified_type Construct_vector_3;

  /*!
  Function object type that provides
  `Triangle_3 operator()(Point_3 a, Point_3 b, Point_3 c)`
  to construct the triangle `(a,b,c)`
  */
  typedef unspecified_type Construct_triangle_3;

  /*!
  Function object type that provides
  `Point_3 operator()(Segment_3 s, int i)`
  to return the source or target of `s`, depending on whether `i` is 0 or 1 respectively, and
  `Point_3 operator()(Triangle_3 t, int i)`
  to return the `i`th vertex of `t`.
  */
  typedef unspecified_type Construct_vertex_3;

  /*!
  Function object type that provides
  `Point_3 operator()(Segment_3 s)`
  to return the source of `s`.
  */
  typedef unspecified_type Construct_source_3;

  /*!
  Function object type that provides
  `Point_3 operator()(Segment_3 s)`
  to return the target of `s`.
  */
  typedef unspecified_type Construct_target_3;

  /*!
  Function object type that provides
  `Point_3 operator()(Point_3 p1, FT w1, Point_3 p2)`
  to compute the barycenter of the points `p1` and `p2` with corresponding weights `w1` and `1-w1`, and
  `Point_3 operator()(Point_3 p1, FT w1, Point_3 p2, FT w2, Point_3 p3, FT w3)`
  to compute the barycenter of the points `p1`, `p2`, and `p3` with corresponding weights `w1`, `w2`, and `w3`.
  */
  typedef unspecified_type Construct_barycenter_3;

  /*!
  Function object type that provides
  `FT operator()(Point_3 p1, Point_3 p2)`
  to compute the squared distance between `p1` and `p2`.
  */
  typedef unspecified_type Compute_squared_distance_3;

  /*!
  Function object type that provides
  `Triangle_2 operator()(Triangle_3 t)`
  which computes a 2-dimensional projection of `t` that preserves edge lengths.
  */
  typedef unspecified_type Construct_triangle_3_to_triangle_2_projection;

  /*!
  Function object type that provides
  `Triangle_2 operator()(Triangle_3 t, std::size_t i, Segment_2 base)`
  which computes a 2-dimensional projection of t that preserves edge lengths,
  such that the `i`th edge of the projection of `t` is incident to `base`.

  \pre The length of the `i`th edge of `t` is equal to the length of `base`
  */
  typedef unspecified_type Construct_triangle_3_along_segment_2_flattening;

  /*!
  Function object type that provides
  `FT operator()(Point_2 x0, Point_2 x1, Point_2 p)`
  which computes the parametric distance of `p` along the
  segment `[x0,x1]`.  That is, it computes `t`, such that
  `p = (1.0 - t)*x0 + t*x1`

  \pre `p` is a point in the segment `[x0,x1]`
  */
  typedef unspecified_type Compute_parametric_distance_along_segment_2;

  /*!
  Function object type that provides
  `Barycentric_coordinate operator()(FT a, FT b, FT c)`
  to introduce a new triangular barycentric coordinate.
  */
  typedef unspecified_type Construct_barycentric_coordinate;

  /*!
  Function object type that provides
  `FT operator(Barycentric_coordinate b, std::size_t i)`
  to get the `i`th weight of barycentric coordinate `b`.
  */
  typedef unspecified_type Construct_barycentric_coordinate_weight;

  /*!
  Function object type that provides
  `Barycentric_coordinate operator()(Triangle_2 t, Point_2 p)`
  which computes the Barycentric location of `p` in `t`.
  */
  typedef unspecified_type Construct_barycentric_coordinate_in_triangle_2;

  /*!
  Function object type that provides
  `Barycentric_coordinate operator()(Triangle_3 t, Point_3 p)`
  which computes the Barycentric location of `p` in `t`.
  */
  typedef unspecified_type Construct_barycentric_coordinate_in_triangle_3;

/// @}

/// \name Predicates
/// @{

  /*!
  Function object type that provides
  `std::pair<CGAL::Surface_mesh_shortest_paths_3::Barycentric_coordinate_type,std::size_t> operator()(Barycentric_coordinate b)`,
  which computes the classification and the associated edge (if applicable) of the coordinate `b`
  \details Returns the pair (`type`, `i`), such that `type` is one of the values of `CGAL::Surface_mesh_shortest_paths_3::Barycentric_coordinate_type`
  - If `type` is `CGAL::Surface_mesh_shortest_paths_3::BARYCENTRIC_COORDINATE_ON_VERTEX`, `i` is the index of that vertex
  - If `type` is `CGAL::Surface_mesh_shortest_paths_3::BARYCENTRIC_COORDINATE_ON_BOUNDARY`, `i` is the index of the non-zero edge
    - 0 if (0,1) are the non-zero coordinates
    - 1 if (1,2) are the non-zero coordinates
    - 2 if (2,0) are the non-zero coordinates
  - Otherwise, the value of `i` is undefined.
  */
  typedef unspecified_type Classify_barycentric_coordinate;

  /*!
  Function object type that provides
  `CGAL::Orientation operator()(Point_2 s1, Point_2 s2, Point_2 s3)`,
  which performs an orientation test for the given points
  */
  typedef unspecified_type Orientation_2;

  /*!
  Function object type that provides
  `CGAL::Comparison_result operator()(Point_2 p1, Point_2 p2, Point_2 p3, Point_2 p4)`, which compares the
  squared distance between `p1` and `p2` with the squared distance between `p3` and `p4`.
  */
  typedef unspecified_type Compare_distance_2;

  /*!
  Function object type that provides
  `CGAL::Comparison_result operator()(Segment_2 s1, Line_2 l1, Segment_2 s2, Line_2 l2)`.
  This compares the relative parametric intersections of `s1` with `l1` against `s2` with `l2`.
  That is, compare the distance of the intersection of `s1` with `l1` from the
  source of `s1`, scaled by the length of `s1`, to the distance of the intersection
  of `s2` with `l2` from the source of `s2`, scaled by the length of `s2`.

  \pre the intersection of `s1` and `l1` is a point
  \pre the intersection of `s2` and `l2` is a point
  */
  typedef unspecified_type Compare_relative_intersection_along_segment_2;

  /*!
  Function object type that provides
  `template <class VertexPointMap> bool operator()(boost::graph_traits<Triangle_mesh>::%vertex_descriptor v, Triangle_mesh& tm, VertexPointMap vpm)`
  that returns true if the vertex is a saddle vertex (more than \f$ 2 \pi \f$ surface area
  over all adjacent faces), and false otherwise.  `vpm` must be a model of concept `ReadablePropertyMap` that maps from `vertex_descriptor` to
  `Point_3` objects.
  */
  typedef unspecified_type Is_saddle_vertex;

/// @}

/// \name Creation
/// @{
  /*!
  */
  SurfaceMeshShortestPathTraits(SurfaceMeshShortestPathTraits& copy);

/// @}

/// \name Operations
/// For all of the above predicate and construction types, e.g. `Func_obj_type`, a function must exist with the name `func_obj_type_object()` that creates an instance of the construction or predicate object type.
/// For example:
/// @{

  Construct_point_2 construct_point_2_object();

/// @}

};