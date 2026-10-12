/*!
\ingroup PkgAlphaWrap2Concepts
\cgalConcept

The concept `AlphaWrapTraits_2` defines the requirements for the geometric traits class
of an alpha wrap oracle.

\cgalRefines{AABBRayIntersectionGeomTraits_2, DelaunayTriangulationTraits_2, PolygonTraits_2}

\cgalHasModelsBegin
\cgalHasModelsBare{Any \cgal %kernel is a model of this traits concept}
\cgalHasModelsEnd
*/

class AlphaWrapTraits_2
{
public:
  /*!
  The field type, must be a model of `FieldNumberType` and `FromDoubleConstructible`
  */
  typedef unspecified_type FT;

  /*!
  A predicate object that must provide the following function operator:

  `bool operator()(Segment_2 s)`,

  which returns `true` iff the segment is degenerate.
  */
  typedef unspecified_type Is_degenerate_2;

  /*!
  A predicate object that must provide the following function operator:

  `Angle operator()(Point_2 p, Point_2 q, Point_2 r)`,

  which returns the angle at point `q` formed by the points `p`, `q`, and `r`.
  */
  typedef unspecified_type Angle_2;

  /*!
  A construction object that must provide the following function operators:

  `CGAL::Bbox_2 operator()(Point_2 p)`,

  `CGAL::Bbox_2 operator()(Segment_2 s)`,

  `CGAL::Bbox_2 operator()(Triangle_2 tr)`,

  that constructs the bounding box of these 2D objects.
  */
  typedef unspecified_type Construct_bbox_2;

  /*!
  A construction object that must provide the following function operator:

  `Vector_2 operator()(Vector_2 v, FT l)`,

  that constructs the scaled vector `l * v`.
  */
  typedef unspecified_type Construct_scaled_vector_2;

  /*!
  A construction object that must provide the following function operator:

  `Point_2 operator()(Point_2 p, Vector_2 v)`,

  that returns the point obtained by translating `p` by the vector `v`.
  */
  typedef unspecified_type Construct_translated_point_2;

  /*!
  A construction object that must be a model of `Kernel::ConstructVertex_2`
  */
  typedef unspecified_type Construct_vertex_2;

  /*!
  A predicate object that must provide the following function operator:

  `bool operator()(Triangle_2 tr, Point_2 p)`,

  which returns whether `p` is on the bounded side of the triangle `tr` or not.
  */
  typedef unspecified_type Has_on_bounded_side_2;

  /*!
  A predicate object that must provide the following function operator:

  `bool operator()(Circle_2 c, Point_2 p)`,

  which returns whether `p` is on the unbounded side of the circle `c` or not.
  */
  typedef unspecified_type Has_on_unbounded_side_2;

  /*!
  A predicate object that must provide the following function operator:

  ` Bounded_side operator()(Point_2 p, Point_2 q, Point_2 t)`,

  which returns the position of the point `t` relative to the circle that has the line segment `pq` as its diameter.
  */
  typedef unspecified_type Side_of_bounded_circle_2;

  // ===

  /*!
  returns the `Is_degenerate_2` predicate.
  */
  Is_degenerate_2 is_degenerate_2_object();

  /*!
  returns the `Angle_2` functor.
  */
  Angle_2 angle_2_object();

  /*!
  returns the `Construct_bbox_2` functor.
  */
  Construct_bbox_2 construct_bbox_2_object();

  /*!
  returns the `Construct_scaled_vector_2` functor.
  */
  Construct_scaled_vector_2 construct_scaled_vector_2_object();

  /*!
  returns the `Construct_translated_point_2` functor.
  */
  Construct_translated_point_2 construct_translated_point_2_object();

  /*!
  returns the `Construct_vertex_2` functor.
  */
  Construct_vertex_2 construct_vertex_2_object();

  /*!
  returns the `Has_on_bounded_side_2` predicate.
  */
  Has_on_bounded_side_2 has_on_bounded_side_2_object();

  /*!
  returns the `Has_on_unbounded_side_2` predicate.
  */
  Has_on_unbounded_side_2 has_on_unbounded_side_2_object();

  /*!
  returns the `Side_of_bounded_circle_2` predicate.
  */
  Side_of_bounded_circle_2 side_of_bounded_circle_2_object();
};
