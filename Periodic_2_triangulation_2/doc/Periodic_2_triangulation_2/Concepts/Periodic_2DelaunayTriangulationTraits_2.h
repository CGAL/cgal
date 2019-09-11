// Copyright (c) 1997-2013 INRIA Sophia-Antipolis (France).
// All rights reserved.

/*!
\ingroup PkgPeriodic2Triangulation2Concepts
\cgalConcept

The concept `Periodic_2DelaunayTriangulationTraits_2` is the first template parameter of the class
`Periodic_2_Delaunay_triangulation_2`. It refines the concepts
`Periodic_2TriangulationTraits_2` and `DelaunayTriangulationTraits_2`.
It redefines the geometric objects, predicates and constructions to work with
point-offset pairs. In most cases the offsets will be (0,0) and the
predicates from `DelaunayTriangulationTraits_2` can be used
directly. For efficiency reasons we maintain for each functor the
version without offsets.

\cgalRefines `DelaunayTriangulationTraits_2` and `Periodic_2TriangulationTraits_2`

In addition to the requirements of the concepts `Periodic_2TriangulationTraits_2`
and `DelaunayTriangulationTraits_2`,
the concept `::Periodic_2DelaunayTriangulationTraits_2` provides a predicate to check the empty circle property. The
corresponding predicate type is called type `::Periodic_2DelaunayTriangulationTraits_2::Side_of_oriented_circle_2`.

The additional constructor object `::Periodic_2DelaunayTriangulationTraits_2::Construct_circumcenter_2` is
used to build the dual Voronoi diagram and are required only if the
dual functions are called. The additional predicate type
`::Periodic_2DelaunayTriangulationTraits_2::Compare_distance_2` is required if calls to
`nearest_vertex(..)` are issued.

\cgalHasModel `CGAL::Periodic_2_Delaunay_triangulation_traits_2<Traits, Offset>`

\sa `DelaunayTriangulationTraits_2`
*/
class Periodic_2DelaunayTriangulationTraits_2
{
public:

/// \name Types
/// @{

  /*!
  Predicate object. Must
  provide the operators

  `Oriented_side operator()(Point p, Point q, Point r, Point s)`
  which takes four points \f$ p, q, r, s\f$ as arguments and returns
  `ON_POSITIVE_SIDE`, `ON_NEGATIVE_SIDE` or,
  `ON_ORIENTED_BOUNDARY` according to the position of points `s`
  with respect to the oriented circle through through \f$ p,q\f$
  and \f$ r\f$ and

  `Oriented_side operator()( Point p, Point q, Point r, Point s, Periodic_2_offset_2 o_p, Periodic_2_offset_2 o_q Periodic_2_offset_2 o_r, Periodic_2_offset_2 o_s)`
  which takes four points \f$ (p, o_p), (q, o_q), (r, o_r), (s, o_s)\f$ as arguments and returns
  `ON_POSITIVE_SIDE`, `ON_NEGATIVE_SIDE` or,
  `ON_ORIENTED_BOUNDARY` according to the position of points `(s, o_s)`
  with respect to the oriented circle through through `(p, o_p), (q, o_q)`
  and `(r, o_r)`.

  This type is required only if the function
  `side_of_oriented_circle(Face_handle f, Point p)` is called.
  */
  typedef unspecified_type Side_of_oriented_circle_2;

  /*!
  Constructor
  object. Provides the operators:

  `Point operator()(Point p, Point q, Point r)`

  which returns
  the circumcenter of the three points `p, q` and `r`.

  `Point operator()(Point p, Point q, Point r, Periodic_2_offset_2 o_p, Periodic_2_offset_2 o_q Periodic_2_offset_2 o_r)`

  which returns the circumcenter of the three points `(p, o_p), (q, o_q)` and `(r, o_r)`.

  This type is required only if the function `Point circumcenter(Face_handle f)`is called.
  */
  typedef unspecified_type Construct_circumcenter_2;

  /*!
  Predicate type. Provides
  the operators:

  `Comparison_result operator()(Point_2 p, Point_2 q, Point_2 r)`
  which returns `SMALLER`, `EQUAL` or `LARGER` according
  to the distance between p and q being smaller, equal or larger than
  the distance between p and r.

  `Comparison_result operator()(Point_2 p, Point_2 q, Point_2 r, Periodic_2_offset_2 o_p, Periodic_2_offset_2 o_q Periodic_2_offset_2 o_r)` which returns `SMALLER`, `EQUAL`
  or `LARGER` according to the distance between `(p, o_p)`,
  and `(q, o_q)` being smaller, equal or larger than
  the distance between `(p, o_p)` and `(r, o_r)`.

  This type is only require if
  `nearest_vertex` queries are issued.
  */
  typedef unspecified_type Compare_distance_2;

/// @}

/// \name Predicate Functions
/// @{

  /*!
  Required only
  if `side_of_oriented_circle` is called
  called.
  */
  Side_of_oriented_circle_2
  side_of_oriented_circle_2_object();

  /*!
  Required only if `circumcenter` is called.
  */
  Construct_circumcenter_2 construct_circumcenter_2_object();

  /*!
  Required only if `compare_distance` is called.
  */
  Compare_distance_2 compare_distance_2_object();

/// @}

}; /* end Periodic_2DelaunayTriangulationTraits_2 */

