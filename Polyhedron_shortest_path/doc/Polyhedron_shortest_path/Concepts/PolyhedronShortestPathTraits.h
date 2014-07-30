/*!
\ingroup PkgPolyhedronShortestPathConcepts

\cgalConcept
 
The concept `PolyhedronShortestPathTraits` describes the required types, 
predicates, and constructions required to execute the 
`Polyhedron_shortest_path` algorithm
 
\cgalHasModel `CGAL::Polyhedron_shortest_path_default_traits<K,P>`
 
*/

class PolyhedronShortestPathTraits
{
public:

/// \name Types
/// @{

  /// The polyhedron type for the algorithm
  typedef unspecified_type Polyhedron;
  
  /// The numeric type for the algorithm.  All geometric types are expected to use this number type.
  /// The square root operation must be supported on this type.
  typedef unspecified_type FT;
  
  /// The 2-dimensional Point type for this algorithm
  typedef unspecified_type Point_2;
  
  /// The 2-dimensional Vector type for this algorithm
  typedef unspecified_type Vector_2;
  
  /// The 2-dimensional Ray type for this algorithm
  typedef unspecified_type Ray_2;
  
  /// The 2-dimensional Line type for this algorithm
  typedef unspecified_type Line_2;

  /// The 2-dimensional Segment type for this algorithm
  typedef unspecified_type Segment_2;
  
  /// The 2-dimensional Triangle type for this algorithm
  typedef unspecified_type Triangle_2;
  
  /// An ordered triple to specify Barycentric coordinates in triangles.
  typedef unspecified_type Barycentric_coordinate;
  
  /// The 3-dimensional Point type for this algorithm
  /// (Must be exactly the same type used by Polyhedron).
  typedef unspecified_type Point_3;
  
  /// The 3-dimensional Vector type for this algorithm
  typedef unspecified_type Vector_3;
  
  /// The 3-dimensional Ray type for this algorithm
  typedef unspecified_type Ray_3;
  
  /// The 3-dimensional Triangle type for this algorithm
  typedef unspecified_type Triangle_3;
  
/// @}

/// \name Constructions
/// @{

  /*!
  A construction object type. Must provide
  `Point_2 operator()(FT x, FT y)` and 
  `Point_2 operator()(CGAL::Origin)` to construct points with
  cartesian coordinates (x,y) and (0,0) respectively.
  */
  typedef unspecified_type Construct_point_2;
  
  /*!
  A construction object type. Must provide
  `Vector_2 operator()(Point_2 a, Point_2 b)`
  to construct the vector b - a
  */
  typedef unspecified_type Construct_vector_2;
  
  /*!
  A construction object type. Must provide
  `Ray_2 operator()(Point_2 a, Point_2 b)`
  to construct the ray originating from a and passing through b
  */
  typedef unspecified_type Construct_ray_2;
  
  /*!
  A construction object type. Must provide
  `Line_2 operator()(Segment_2 s)` and
  `Line_2 operator()(Ray_2 r)`
  to construct the supporting line to s and r respectively.
  */
  typedef unspecified_type Construct_line_2;
  
  /*!
  A construction object type. Must provide
  `Segment_2 operator()(Point_2 a, Point_2 b)`
  to construct the segment (a,b)
  */
  typedef unspecified_type Construct_segment_2;
  
  /*!
  A construction object type. Must provide
  `Triangle_2 operator()(Point_2 a, Point_2 b, Point_2 c)`
  to construct the triangle (a,b,c)
  */
  typedef unspecified_type Construct_triangle_2;

  /*!
  A construction object type. Must provide
  `Point_2 operator()(Segment_2 s, int i)`
  to return the source or target of s, depending on whether i is 0 or 1 respectively, and
  `Point_2 operator()(Triangle_2 t, int i)`
  to return the ith vertex of t.
  */
  typedef unspecified_type Construct_vertex_2;
  
  /*!
  A construction object type. Must provide
  `Point_2 operator()(Segment_2 s)`
  to return the source of s.
  */
  typedef unspecified_type Construct_source_2;
  
  /*!
  A construction object type. Must provide
  `Point_2 operator()(Segment_2 s)`
  to return the target of s.
  */
  typedef unspecified_type Construct_target_2;
  
  /*!
  A construction object type. Must provide
  `Point_2 operator()(Point_2 p1, FT w1, Point_2 p2, FT w2)`
  to compute the barycenter of the points p1 and p2 with corresponding weights w1 and w2, and
  `Point_2 operator()(Point_2 p1, FT w1, Point_2 p2, FT w2, Point_2 p3, FT w3)`
  to compute the barycenter of the points p1, p2, and p3 with corresponding weights w1, w2, and w3.
  */
  typedef unspecified_type Construct_barycenter_2;
  
  /*!
  A construction object type. Must provide
  `FT operator()(A obj1, B obj2)`
  to compute the squared distance between obj1 and obj2, A and B can be any one of 
  Point_2,
  Segment_2,
  type objects.
  */
  typedef unspecified_type Compute_squared_distance_2;
  
  /*!
  A construction object type. 
  Must support the result_of protocol, that is the return type of the operator()(A, B) is CGAL::cpp11::result<Intersect_2(A,B)>.
  Must provide
  `CGAL::cpp11::result<Intersect_2(A,B)> operator()(A obj1, B obj2)`
  to compute the intersection between obj1 and obj2, where A and B can be any one of
  Line_2,
  Ray_2,
  Segment_2,
  type objects.
  */
  typedef unspecified_type Intersect_2;
  
  /*!
  A construction object type. Must provide
  `Point_2 operator()(Ray_2 r, int i)`
  to construct a point on r, such that i = 0 is the source, and any i > 0 is not the source.
  */
  typedef unspecified_type Construct_point_on_2;

  /*!
  A construction object type. Must provide
  `Vector_3 operator()(Point_3 a, Point_3 b)`
  to construct the vector b - a
  */
  typedef unspecified_type Construct_vector_3;

  /*!
  A construction object type. Must provide
  `Triangle_3 operator()(Point_3 a, Point_3 b, Point_3 c)`
  to construct the triangle (a,b,c)
  */
  typedef unspecified_type Construct_triangle_3;
  
  /*!
  A construction object type. Must provide
  `Point_3 operator()(Segment_3 s, int i)`
  to return the source or target of s, depending on whether i is 0 or 1 respectively, and
  `Point_3 operator()(Triangle_3 t, int i)`
  to return the ith vertex of t.
  */
  typedef unspecified_type Construct_vertex_3;
  
  /*!
  A construction object type. Must provide
  `Point_3 operator()(Segment_3 s)`
  to return the source of s.
  */
  typedef unspecified_type Construct_source_3;
  
  /*!
  A construction object type. Must provide
  `Point_3 operator()(Segment_3 s)`
  to return the target of s.
  */
  typedef unspecified_type Construct_target_3;
  
  /*!
  A construction object type. Must provide
  `Point_3 operator()(Point_3 p1, FT w1, Point_3 p2)`
  to compute the barycenter of the points p1 and p2 with corresponding weights w1 and 1-w1, and
  `Point_3 operator()(Point_3 p1, FT w1, Point_3 p2, FT w2, Point_3 p3, FT w3)`
  to compute the barycenter of the points p1, p2, and p3 with corresponding weights w1, w2, and w3.
  */
  typedef unspecified_type Construct_barycenter_3;
  
  /*!
  A construction object type. Must provide
  `FT operator()(Point_3 p1, Point_3 p2)`
  to compute the squared distance between p1 and p2.
  */
  typedef unspecified_type Compute_squared_distance_3;

  /*!
  A construction object type.  Must provide
  `Triangle_2 operator()(Triangle_3 t)`
  which computes a 2-dimensional projection of t
  that preserves all edge lengths.
  */
  typedef unspecified_type Project_triangle_3_to_triangle_2;
  
  /*!
  A construction object type.  Must provide
  `Triangle_2 operator()(Triangle_3 t, size_t i, Segment_2 base)`
  which computes a 2-dimensional projection of t
  that preserves all edge lengths, such that the `i`th edge of t 
  lies along `base`.
  
  \pre The length of the `i`th edge of t is equal to the length of `base`
  */
  typedef unspecified_type Flatten_triangle_3_along_segment_2;
  
  /*! 
  A construction object type.  Must provide
  `FT operator()(Point_2 x0, Point_2 x1, Point_2 p)`
  which computes the parametric distance of p along the 
  segment [x0,x1].  That is, it computes `t`, such that
  p = (1.0 - t)*x0 + t*x1
  
  \pre p is a point along the segment [x0,x1]
  */
  typedef unspecified_type Parametric_distance_along_segment_2;

  /*!
  A construction object type.  Must provide
  `Barycentric_coordinate operator()(FT a, FT b, FT c)'
  to introduce a new triangular barycentric coordinate.
  */
  typedef unspecified_type Construct_barycentric_coordinate;
  
  /*!
  A construction object type.  Must provide
  `FT operator(Barycentric_coordinate b, size_t i)'
  to get the `i`th weight of barycentric coordinate `b`.
  */
  typedef unspecified_type Construct_barycentric_coordinate_weight;
  
  /*!
  A construction object type.  Must provide
  `Barycentric_coordinate operator()(Triangle_2 t, Point_2 p)`
  which computes the Barycentric location of p in t.
  */
  typedef unspecified_type Construct_barycentric_coordinate_in_triangle_2;
  
  /*!
  A construction object type.  Must provide
  `Barycentric_coordinate operator()(Triangle_3 t, Point_3 p)`
  which computes the Barycentric location of p in t.
  */
  typedef unspecified_type Construct_barycentric_coordinate_in_triangle_3;
  
/// @}

/// \name Predicates
/// @{

  /*!
  A predicate object type.  Must provide:
  `std::pair<CGAL::internal::Barycentric_coordinate_type,size_t> operator()(Barycentric_coordinate b)`,
  which computes the classification and nearest edge of the coordinate `b`
  \details Returns the pair (`type`, `i`), such that `type` is:
  - `CGAL::internal::BARYCENTRIC_COORDINATE_VERTEX` if `b` has exactly one non-zero component equal to 1, and the rest are zero
  - `CGAL::internal::BARYCENTRIC_COORDINATE_EDGE` if `b` has exactly one zero component, and the rest sum to 1
  - `CGAL::internal::BARYCENTRIC_COORDINATE_INTERNAL` if `b` has no non-zero component, and they all sum to 1
  - `CGAL::internal::BARYCENTRIC_COORDINATE_EXTERNAL` if the components of `b` do not sum to 1
  - `CGAL::internal::BARYCENTRIC_COORDINATE_INVALID` otherwise
  If `type` is `CGAL::internal::BARYCENTRIC_COORDINATE_VERTEX`, `i` should be the index of that vertex
  If `type` is `CGAL::internal::BARYCENTRIC_COORDINATE_EDGE`, `i` should be the index of the non-zero edge
  - 0 if (0,1) are the non-zero coordinates
  - 1 if (1,2) are the non-zero coordinates
  - 2 if (2,0) are the non-zero coordinates
  otherwise, the value of `i` does not matter.
  */
  typedef unspecified_type Classify_barycentric_coordinate;
  
  /*!
  A predicate object type. Must provide 
  `CGAL::Orientation operator()(Point_2 s1, Point_2 s2, Point_2 s3)`,
  which performs an orientation test for the given points
  */
  typedef unspecified_type Orientation_2;
  
  /*!
  A predicate object type. Must provide
  `CGAL::Comparison_result operator()(Point_2 p1, Point_2 p2, Point_2 p3, Point_2 p4)`, which compares the
  squared distance between `p1` and `p2` with the squared distance between `p3` and `p4`.
  */
  typedef unspecified_type Compare_distance_2;
  
  /*!
  A predicate object type.  Must provide
  `CGAL::Comparison_result operator()(Segment_2 s1, Line_2 l1, Segment_2 s2, Line_2 l2)`.
  This compares the relative parametric intersections of s1 with l1 against s2 with l2.
  That is, compare the distance of the intersection of s1 with l1 from the 
  start point of s1, scaled by the length of s1, to the distance of the intersection 
  of s2 with l2 from the start point of s2, scaled by the length of s2.
  
  \pre the intersection of s1 and l1 is a single point
  \pre the intersection of s2 and l2 is a single point
  */
  typedef unspecified_type Compare_relative_intersection_along_segment_2;
  
  /*!
  A predicate object type.  Must provide 
  `template <class VertexPointMap> bool operator()(boost::graph_traits<Polyhedron>::vertex_descriptor v, Polyhedron& p, VertexPointMap vpm)`
  that returns true if the vertex is a saddle vertex (more than 2pi surface area 
  over all adjacent faces), and false otherwise.  vmp must be a a model of concept ReadablePropertyMap that maps from vertex_descriptor to
  Point_3 objects.
  */
  typedef unspecified_type Is_saddle_vertex;

/// @}

/// \name Methods
/// @{

  /*!
  \brief returns an instance of Classify_barycentric_coordinate
  */
  Classify_barycentric_coordinate classify_barycentric_coordinate_object();
  
  /*!
  \brief returns an instance of Orientation_2
  */
  Orientation_2 orientation_2_object();
  
  /*!
  \brief returns an instance of Intersect_2
  */
  Intersect_2 intersect_2_object();
  
  /*!
  \brief returns an instance of Construct_point_2
  */
  Construct_point_2 construct_point_2_object();
  
  /*!
  \brief returns an instance of Construct_vector_2
  */
  Construct_vector_2 construct_vector_2_object();

  /*!
  \brief returns an instance of Construct_ray_2
  */
  Construct_ray_2 construct_ray_2_object();
  
  /*!
  \brief return an instance of Construct_line_2
  */
  Construct_line_2 construct_line_2_object();
  
  /*!
  \brief return an instance of Construct_segment_2
  */
  Construct_segment_2 construct_segment_2_object();
  
  /*!
  \brief return an instance of Construct_triangle_2
  */
  Construct_triangle_2 construct_triangle_2_object();
  
  /*!
  \brief returns an instance of Compute_squared_distance_2
  */
  Compute_squared_distance_2 compute_squared_distance_2_object();
  
  /*!
  \brief returns an instance of Construct_vector_3
  */
  Construct_vector_3 construct_vector_3_object();
  
  /*!
  \brief returns an instance of Construct_triangle_3
  */
  Construct_triangle_3 construct_triangle_3_object();
  
  /*!
  \brief returns an instance of Construct_triangle_3
  */
  Construct_vertex_3 construct_vertex_3_object();
  
  /*!
  \brief returns an instance of Construct_source_3
  */
  Construct_source_3 construct_source_3_object();
  
  /*!
  \brief returns an instance of Construct_target_3
  */
  Construct_target_3 construct_target_3_object();
  
  /*!
  \brief returns an instance of Construct_barycenter_3
  */
  Construct_barycenter_3 construct_barycenter_3_object();

  /*!
  \brief returns an instance of Compute_squared_distance_3
  */
  Compute_squared_distance_3 compute_squared_distance_3_object();
  
  /*!
  \brief returns an instance of Compare_relative_intersection_along_segment_2
  */
  Compare_relative_intersection_along_segment_2 compare_relative_intersection_along_segment_2_object();
  
  /*!
  \brief returns an instance of Construct_projected_point_2
  */
  Is_saddle_vertex is_saddle_vertex_object();
  
  /*!
  \brief returns an instance of Project_triangle_3_to_triangle_2
  */
  Project_triangle_3_to_triangle_2 project_triangle_3_to_triangle_2_object();
  
  /*!
  \brief returns an instance of Flatten_triangle_3_along_segment_2
  */
  Flatten_triangle_3_along_segment_2 flatten_triangle_3_along_segment_2_object();
  
  /*!
  \brief returns an instance of Construct_barycentric_coordinate
  */
  Construct_barycentric_coordinate construct_barycentric_coordinate_object();
  
  /*!
  \brief returns an instance of Construct_barycentric_coordinate_weight
  */
  Construct_barycentric_coordinate_weight construct_barycentric_coordinate_weight_object();
  
  /*!
  \brief returns an instance of Construct_barycentric_coordinate_in_triangle_2
  */
  Construct_barycentric_coordinate_in_triangle_2 construct_barycentric_coordinate_in_triangle_2_object();

  /*!
  \brief returns an instance of Construct_barycentric_coordinate_in_triangle_3
  */
  Construct_barycentric_coordinate_in_triangle_3 construct_barycentric_coordinate_in_triangle_3_object();
  
  /*!
  \brief returns an instance of Parametric_distance_along_segment_2
  */
  Parametric_distance_along_segment_2 parametric_distance_along_segment_2_object();

/// @}

};