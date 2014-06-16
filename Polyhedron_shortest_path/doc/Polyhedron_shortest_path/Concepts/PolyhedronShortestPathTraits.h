/*!
\ingroup PkgPolyhedronShortestPathConcepts

\cgalConcept
 
The concept `PolyhedronShortestPathTraits' describes the required types, 
predicates, and constructions required to execute the 
'Polyhedron_shortest_path' algorithm
 
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
  
  /// An ordered triple to specify Barycentric in triangles.
  /// Must support member-wise addition, and scalar multiplication.
  /// (therefore, and alias of Vector_3 is sufficient)
  typedef unspecified_type Barycentric_coordinate;
  
  /// The 3-dimensional Point type for this algorithm
  /// (Must be exactly the same type used by Polyhedron).
  typedef unspecified_type Point_3;
  
  /// The 3-dimensional Vector type for this algorithm
  typedef unspecified_type Vector_3;
  
  /// The 3-dimensional Triangle type for this algorithm
  typedef unspecified_type Triangle_3;
  
/// @}

/// \name Predicates
/// @{

  /*!
  A predicate object type. Must provide 
  `CGAL::Orientation operator()(Point_2 s1, Point_2 s2, Point_2 s3)`,
  which performs an orientation test for the given points
  */
  typedef unspecified_type Orientation_2;
  
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

/// \name Constructions
/// @{

  /*!
  A construction object type.  Must provide an intersection 
  method between all possible pairs of 2D primitives required
  by this concept
  */
  typedef unspecified_type Intersection_2;
  
  /*!
  A construction object type.  Must provide 
  `FT operator()(Point_2 p1, Point_2 p2)`
  which computes the squared euclidean distance between p1 and p2.
  */
  typedef unspecified_type Compute_squared_distance_2;
  
  /*!
  A construction object type.  Must provide 
  `FT operator()(Point_3 p1, Point_3 p2)`
  which computes the squared euclidean distance between p1 and p2.
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
  typedef unspecified_type Parameteric_distance_along_segment_2;
  
  /*!
  A construction object type.  Must provide
  `Barycentric_coordinate operator()(Triangle_2 t, Point_2 p)`
  which computes the Barycentric location of p in t.
  */
  typedef unspecified_type Construct_barycentric_coordinate_2;
  
  /*!
  A construction object type.  Must provide
  `Barycentric_coordinate operator()(Triangle_3 t, Point_3 p)`
  which computes the Barycentric location of p in t.
  */
  typedef unspecified_type Construct_barycentric_coordinate_3;
  
  /*!
  A construction object type.  Must provide
  `Point_2 operator()(Triangle_2 t, Barycentric_coordinate a)`
  which computes the location of the barycentric coordinate a
  in triangle t.  
  */
  typedef unspecified_type Construct_triangle_location_2;
  
  /*!
  A construction object type.  Must provide
  `Point_3 operator()(Triangle_3 t, Barycentric_coordinate a)`
  which computes the location of the barycentric coordinate a
  in triangle t.  
  */
  typedef unspecified_type Construct_triangle_location_3;
  
/// @}

/// \name Methods
/// @{

  /*!
  \brief returns an instance of the Orientation_2 object
  */
  Orientation_2 orientation_2_object();
  
  /*!
  \brief returns an instance of Intersect_2
  */
  Intersect_2 intersect_2_object();
  
  /*!
  \brief returns an instance of Construct_projected_point_2
  */
  Construct_projected_point_2 construct_projected_point_2_object();
  
  /*!
  \brief returns an instance of Compute_squared_distance_2
  */
  Compute_squared_distance_2 compute_squared_distance_2_object();
  
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
  \brief returns an instance of Construct_barycentric_coordinate_2
  */
  Construct_barycentric_coordinate_2 construct_barycentric_coordinate_2_object();
  
  /*!
  \brief returns an instance of Construct_triangle_location_2
  */
  Construct_triangle_location_2 construct_triangle_location_2_object();
  
  /*!
  \brief returns an instance of Construct_barycentric_coordinate_3
  */
  Construct_barycentric_coordinate_3 construct_barycentric_coordinate_3_object();
  
  /*!
  \brief returns an instance of Construct_triangle_location_3
  */
  Construct_triangle_location_3 construct_triangle_location_3_object();
  
  /*!
  \brief returns an instance of Parameteric_distance_along_segment_2
  */
  Parameteric_distance_along_segment_2 parameteric_distance_along_segment_2_object();

/// @}

};