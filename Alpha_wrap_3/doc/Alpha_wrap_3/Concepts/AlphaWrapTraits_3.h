// Undocumented for now
#ifndef DOXYGEN_RUNNING

/*!
\ingroup PkgAlphaWrap3Concepts
\cgalConcept

@fixme Because of a few calls to PMP::, if you only look at the doc, it's all pointless because
you require Kernel. Stitch_borders doesn't even have clear geometric traits requirements...

The concept `AlphaWrapTraits_3` defines the requirements for the geometric traits class
of an alpha wrap oracle.

\cgalHasModelsBegin
\cgalHasModelsBare{Any 3D %kernel is a model of this traits concept}
\cgalHasModelsEnd
*/

class AlphaWrapTraits_3
{
public:
  /// The field type
  typedef unspecified_type FT;

  /// The point type
  typedef unspecified_type Point_3;

  /// The vector type
  typedef unspecified_type Vector_3;

  /// The triangle type
  typedef unspecified_type Triangle_3;

  /// The tetrahedron type
  typedef unspecified_type Tetrahedron_3;

  /// The ball type
  typedef unspecified_type Ball_3;

  /*!
  A predicate object that must provide the following function operator:

  `CGAL::Comparison_result operator()(Point_3 p, Point_3 q, FT sqd)`,

  which compares the squared distance between the two points `p` and `q` to the value `sqd`.
  */
  typedef unspecified_type Compare_squared_distance_3;

  /*!
  A construction object that must provide the following function operator:

  `FT operator()(Vector_3 v)`,

  which returns the squared length of the vector `v`.
  */
  typedef unspecified_type Compute_squared_length_3;

  /*!
  A construction object that must provide the following function operators:

  `FT operator()(Point_3 p, Point_3 q, Point_3 r)`,

  and

  `FT operator()(Point_3 p, Point_3 q, Point_3 r, Point_3 s)`,

  which return the squared radius of the smallest sphere enclosing the points.
  */
  typedef unspecified_type Compute_squared_radius_3;

  /*!
  A predicate object that must provide the following function operator:

  `bool operator()(Point_3 p, Point_3 q, Point_3 r)`,

  which returns `true`, iff `p`, `q`, and `r` are collinear.
  */
  typedef unspecified_type Collinear_3;

  /*!
  A construction object that must provide the following function operator:

  `Ball_3 operator()(Point_3 p, FT sqr)`,

  which returns the ball centered at `p` with squared radius `sqr`.
  */
  typedef unspecified_type Construct_ball_3;

  /*!
  A construction object that must provide the following function operators:

  `CGAL::Bbox_3 operator()(Triangle_3 tr)`,

  and

  `CGAL::Bbox_3 operator()(Tetrahedron_3 tet)`,

  which return an axis-aligned bounding box that encloses the object.
  */
  typedef unspecified_type Construct_bbox_3;

  /*!
  A construction object that must provide the following function operator:

  `Vector_3 operator()(Vector_3 v, FT a)`,

  which returns the vector `v` with length scaled by `a`.
  */
  typedef unspecified_type Construct_scaled_vector_3;

  /*!
  A construction object that must provide the following function operator:

  `Point_3 operator()(Point_3 p, Vector_3 v)`,

  which returns the point `p` translated by the vector `v`.
  */
  typedef unspecified_type Construct_translated_point_3;

  /*!
  A predicate object that must provide the following function operators:

  `bool operator()(Tetrahedron_3 tet, Point_3 p)`,

  and

  `bool operator()(Sphere_3 s, Point_3 p)`,

  which return `true` iff `p` is on the bounded side of the respective kernel objects.
  */
  typedef unspecified_type Has_on_bounded_side_3;

  /*!
  A predicate object that must provide the following function operator:

  `bool operator()(Point_3 p, Point_3 q, Point_3 r, Point_3 s)`,

  which returns `true` iff `s` is on the bounded side of the smallest sphere enclosing `p`, `q`, and `r`.
  */
  typedef unspecified_type Side_of_bounded_sphere_3;

  // ===

  /*!
  returns the `Compare_squared_distance_3` predicate.
  */
  Compare_squared_distance_3 compare_squared_distance_3_object();

  /*!
  returns the `Compute_squared_length_3` predicate.
  */
  Compute_squared_length_3 compute_squared_length_3_object();

  /*!
  returns the `Compute_squared_radius_3` predicate.
  */
  Compute_squared_radius_3 compute_squared_radius_3_object();

  /*!
  returns the `Collinear_3` predicate.
  */
  Collinear_3 collinear_3_object();

  /*!
  returns the `Construct_ball_3` construction.
  */
  Construct_ball_3 construct_ball_3_object();

  /*!
  returns the `Construct_scaled_vector_3` construction.
  */
  Construct_scaled_vector_3 construct_scaled_vector_3_object();

  /*!
  returns the `Construct_translated_point_3` construction.
  */
  Construct_translated_point_3 construct_translated_point_3_object();

  /*!
  returns the `Has_on_bounded_side_3` predicate.
  */
  Has_on_bounded_side_3 has_on_bounded_side_3_object();

  /*!
  returns the `Side_of_bounded_sphere_3` predicate.
  */
  Side_of_bounded_sphere_3 side_of_bounded_sphere_3_object();
};

#endif // DOXYGEN_RUNNING
