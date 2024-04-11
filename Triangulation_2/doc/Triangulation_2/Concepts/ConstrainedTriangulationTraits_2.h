
/*!
\ingroup PkgTriangulation2Concepts
\cgalConcept

The concept `ConstrainedTriangulationTraits_2` defines the requirements for the geometric
traits class of a constrained triangulation
( `CGAL::Constrained_triangulation_2<Traits,Tds,Itag>`)
that supports intersections of input constraints (i. e.
when the template parameter `Itag` is instantiated
by one of the tag classes `CGAL::Exact_intersections_tag` or
`CGAL::Exact_predicates_tag`). This concept refines the concept
`TriangulationTraits_2`, adding requirements for function objects
to compute the intersection points of two constraints.
When `CGAL::Exact_predicates_tag` is used, the
traits class is
also required to provide additional types
to compute the squared distance between a point and a line

\cgalRefines{TriangulationTraits_2}

\cgalHasModelsBegin
\cgalHasModelsBare{All models of the \cgal concept `Kernel`}
\cgalHasModels{CGAL::Projection_traits_3<K>}
\cgalHasModels{CGAL::Projection_traits_xy_3<K>}
\cgalHasModels{CGAL::Projection_traits_yz_3<K>}
\cgalHasModels{CGAL::Projection_traits_xz_3<K>}
\cgalHasModelsEnd

\sa `TriangulationTraits_2`
\sa `ConstrainedDelaunayTriangulationTraits_2`
\sa `CGAL::Constrained_triangulation_2<Traits,Tds,Itag>`

*/

class ConstrainedTriangulationTraits_2 {
public:

/// \name Types
/// @{

/*!
A function object whose `operator()` computes the intersection of two segments.

`std::optional<std::variant<Point_2,Segment_2> > operator()(Segment_2 s1, Segment_2 s2);`
Returns the intersection of `s1` and `s2`.
*/
typedef unspecified_type Intersect_2;

///@}

/// \name Types required with Exact_predicates_tag
/// When the constrained triangulation is instantiated with the intersection tag `CGAL::Exact_predicates_tag`, the used algorithm needs to be able to compare some distances between points and lines and the following types are further required.
/// @{

/*!
A number type supporting the comparison operator `<`.
*/
typedef unspecified_type RT;

/*!
The line type.
*/
typedef unspecified_type Line_2;

/*!
A function object whose `operator()`
constructs a line from two points.

`Line_2 operator()(Point_2 p1, Point_2 p2)`.
*/
typedef unspecified_type Construct_line_2;

/*!
A function object whose
`operator()` computes the squared distance between
a line and a point.

`RT operator()(Line_2 l, Point_2 p);` Returns the squared distance
between `p` and `l`.
*/
typedef unspecified_type Compute_squared_distance_2;

/*!
A function object whose
`operator()` computes the bounding box of a point.

`CGAL::Bbox_2 operator()(Point_2 p);` Returns the bounding box of `p`.
The result type is `CGAL::Bbox_2` (even for projection traits classes).
*/
typedef unspecified_type Compute_bounding_box_2;

/// @}

/// \name Access to Constructor Objects
/// @{

/*!

*/
Intersect_2 intersect_2_object();

/*!
required when
the intersection tag is `CGAL::Exact_predicates_tag`.
*/
Construct_line_2 construct_line_2_object();

/*!
required when
the intersection tag is `CGAL::Exact_predicates_tag`.
*/
Compute_squared_distance_2
compute_squared_distance_2_object();

/// @}

}; /* end ConstrainedTriangulationTraits_2 */
