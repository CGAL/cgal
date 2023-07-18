// Copyright (c) 1997-2013 INRIA Sophia-Antipolis (France).
// All rights reserved.

/*!
\ingroup PkgPeriodic2Triangulation2Concepts
\cgalConcept

The concept `Periodic_2TriangulationTraits_2` is the first template parameter of the classes
`Periodic_2_triangulation_2<Traits, Tds>`. This concept provides the types of
the geometric primitives used in the triangulation and some function
object types for the required predicates on those primitives.

It refines the concept
`TriangulationTraits_2` from the \cgal \ref PkgTriangulation2Ref package. It redefines the
geometric objects, predicates and constructions to work with
point-offset pairs. In most cases the offsets will be (0,0) and the
predicates from `TriangulationTraits_2` can be used
directly. For efficiency reasons we maintain for each functor the
version without offsets.

\cgalRefines{TriangulationTraits_2}
In addition to the requirements described for the traits class
`TriangulationTraits_2`, the geometric traits class of a
Periodic triangulation must fulfill the following
requirements:

\cgalHasModelsBegin
\cgalHasModels{CGAL::Periodic_2_triangulation_traits_2}
\cgalHasModelsEnd

\sa `TriangulationTraits_2`
\sa `CGAL::Periodic_2_triangulation_2<Traits,Tds>`

*/

class Periodic_2TriangulationTraits_2
{
public:

/// \name Types
/// @{

  /*!
  The point type. It must be a model of
  `Kernel::Point_2`.
  */
  typedef unspecified_type Point_2;

  /*!
  The segment type. It must be a model
  of `Kernel::Segment_2`.
  */
  typedef unspecified_type Segment_2;

  /*!
  The vector type. It must be a model
  of `Kernel::Vector_2`.
  */
  typedef unspecified_type Vector_2;

  /*!
  The triangle type. It must be a
  model of `Kernel::Triangle_2`.
  */
  typedef unspecified_type Triangle_2;

  /*!
  A type representing an
  axis-aligned rectangle. It must be a model of
  `Kernel::Iso_rectangle_2`.
  */
  typedef unspecified_type Iso_rectangle_2;

  /*!
  The offset type. It must
  be a model of the concept `Periodic_2Offset_2`.
  */
  typedef unspecified_type Periodic_2_offset_2;

/// @}

/// \name Predicate types
/// @{

  /*!

  A predicate object that must provide the function operators

  `Comparison_result operator()(Point_2 p, Point_2 q)`,

  which returns `EQUAL` if the \f$ x\f$-coordinates of the two points are equal and

  `Comparison_result operator()(Point_2 p, Point_2 q, Periodic_2_offset_2 o_p, Periodic_2_offset_2 o_q)`,

  which returns `EQUAL` if the \f$ x\f$-coordinates and \f$ x\f$-offsets of
  the two point-offset pairs are equal. Otherwise it must return a
  consistent order for any two points. \pre `p`, `q` lie inside the domain.
  */
  typedef unspecified_type Compare_x_2;

  /*!

  A predicate object that must provide the function operators

  `Comparison_result operator()(Point_2 p, Point_2 q)`,

  which returns `EQUAL` if the \f$ y\f$-coordinates of the two points are equal and

  `Comparison_result operator()(Point_2 p, Point_2 q, Periodic_2_offset_2 o_p, Periodic_2_offset_2 o_q)`,

  which returns `EQUAL` if the \f$ y\f$-coordinates and \f$ y\f$-offsets of
  the two point-offset pairs are equal. Otherwise it must return a
  consistent order for any two points. \pre `p`, `q` lie inside the domain.
  */
  typedef unspecified_type Compare_y_2;

  /*!

  Predicate object. Provides the operators:

  `bool operator()(Point p, Point q)` and

  `bool operator()(Point p, Point q, Periodic_2_offset_2 o_p, Periodic_2_offset_2 o_q)`

  which returns `true` if `p` is before `q`
  according to the \f$ x\f$-ordering of points.

  This predicate is only necessary if the insert function with a range
  of points (using Hilbert sorting) is used.
  */
  typedef unspecified_type Less_x_2;

  /*!
  Predicate object. Provides
  the operators:

  `bool operator()(Point p, Point q)` and

  `bool operator()(Point p, Point q, Periodic_2_offset_2 o_p, Periodic_2_offset_2 o_q)`

  which returns `true` if `p` is before `q`
  according to the \f$ y\f$-ordering of points.

  This predicate is only necessary if the insert function with a range of
  points (using Hilbert sorting) is used.

  */
  typedef unspecified_type Less_y_2;

  /*!
  A predicate object that must provide the function operators

  `Orientation operator()(Point_2 p, Point_2 q, Point_2 r)`,

  which returns `LEFT_TURN`, `RIGHT_TURN` or `COLLINEAR`
  depending on \f$ r\f$ being, with respect to the oriented line `pq`,
  on the left side, on the right side or on the line.
  and

  `Orientation operator()(Point_2 p, Point_2 q, Point_2 r, Periodic_2_offset_2 o_p, Periodic_2_offset_2 o_q, Periodic_2_offset_2 o_r)`,

  which returns `LEFT_TURN`, `RIGHT_TURN` or `COLLINEAR`
  depending on `(r,o_r)` being, with respect to the oriented line
  defined by `(p,o_p)(q,o_q)` on the left side, on the right side
  or on the line.
  */
  typedef unspecified_type Orientation_2;

/// @}

/// \name Constructor types:
/// Note that the traits must provide exact constructions in order to
/// guarantee exactness of the following construction functors.
/// @{

  /*!
  A constructor object for
  `Point_2`. Provides:

  `Point_2 operator()(Point_2 p,Periodic_2_offset_2 p_o)`,

  which constructs a point from a point-offset pair.
  \pre `p` lies inside the domain.

  */
  typedef unspecified_type Construct_point_2;

  /*!
  A constructor object for
  `Segment_2`. Provides:

  `Segment_2 operator()(Point_2 p,Point_2 q)`,

  which constructs a segment from two points and

  `Segment_2 operator()(Point_2 p,Point_2 q, Periodic_2_offset_2 o_p, Periodic_2_offset_2 o_q)`,

  which constructs a segment from the points `(p,o_p)` and `(q,o_q)`.

  */
  typedef unspecified_type Construct_segment_2;

  /*!
  A constructor object for
  `Triangle_2`. Provides:

  `Triangle_2 operator()(Point_2 p,Point_2 q,Point_2 r )`,

  which constructs a triangle from three points and

  `Triangle_2 operator()(Point_2 p,Point_2 q,Point_2 r, Periodic_2_offset_2 o_p, Periodic_2_offset_2 o_q, Periodic_2_offset_2 o_r)`,

  which constructs a triangle from the three points `(p,o_p)`,
  `(q,o_q)` and `(r,o_r)`.
  */
  typedef unspecified_type Construct_triangle_2;

/// @}

/// \name Creation
/// @{

  /*!
  Default constructor.
  */
  Periodic_2_triangulation_traits_2();

  /*!
  Copy constructor.
  */
  Periodic_2_triangulation_traits_2(const Periodic_2_triangulation_traits_2& tr);

/// @}

/// \name Predicate functions
/// The following functions give access to the predicate and constructor objects.
/// @{

  /*!

  */
  Compare_x_2 compare_x_2_object();

  /*!

  */
  Compare_y_2 compare_y_2_object();

  /*!

  */
  Less_x_2 less_x_2_object();

  /*!

  */
  Less_y_2 less_y_2_object();

  /*!

  */
  Orientation_2 orientation_2_object();

  /*!

  */
  Construct_point_2 construct_point_2_object();

  /*!

  */
  Construct_segment_2 construct_segment_2_object();

  /*!

  */
  Construct_triangle_2 construct_triangle_2_object();

/// @}

/// \name Access Functions
/// @{

  /*!
  Sets the fundamental domain. This is necessary to evaluate predicates
  correctly.
  \pre `domain` represents a square.
  */
  void set_domain(Iso_rectangle_2 domain);

  /*!
  Returns the fundamental domain.
  */
  Iso_rectangle_2 get_domain() const;

/// @}

}; /* end Periodic_2TriangulationTraits_2 */

