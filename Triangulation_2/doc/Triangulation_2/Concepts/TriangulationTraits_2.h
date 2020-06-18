
/*!
\ingroup PkgTriangulation2Concepts
\cgalConcept

\cgalRefines SpatialSortingTraits_2

The concept `TriangulationTraits_2` describes the set of requirements
to be fulfilled by any class used to instantiate the first template
parameter of the class `CGAL::Triangulation_2<Traits,Tds>`.  This concept
provides the types of the geometric primitives used in the
triangulation and some function object types for the required
predicates on those primitives.

\cgalHasModel All models of `Kernel`.
\cgalHasModel `CGAL::Projection_traits_xy_3<K>`
\cgalHasModel `CGAL::Projection_traits_yz_3<K>`
\cgalHasModel `CGAL::Projection_traits_xz_3<K>`

\sa `CGAL::Triangulation_2`

*/

class TriangulationTraits_2 {
public:

/// \name Types
/// @{

/*!
The point type.
*/
typedef unspecified_type Point_2;

/*!
The segment type.
*/
typedef unspecified_type Segment_2;

/*!
The triangle type.
*/
typedef unspecified_type Triangle_2;

/*!
A function object to construct a `Point_2`.

Provides:

`Point_2 operator()(Point_2 p)`,

which simply returns p.

\note It is advised to return a const reference to `p` to avoid useless copies.

\note This peculiar requirement is necessary because `CGAL::Triangulation_2`
internally manipulates points with a `Point` type that is not always `Point_2`.
*/
/*
For example, `CGAL::Regular_triangulation_2` inherits `CGAL::Triangulation_2`
with `Point` being a two-dimensional weighted point. Since some predicates and
constructors (such as `Orientation_2`) can only use `Point_2` objects in arguments,
it is necessary to convert objects of type `Point` to objects of type `Point_2`
before calling these functions, using the kernel functor `Construct_point_2`.
In the setting of a basic triangulation, `Point` and `Point_2` are identical and
so `Construct_point_2` is simply the identity. Refinements of this concept will
require more significant overloads to the `Construct_point_2` functor.
*/
typedef unspecified_type Construct_point_2;

/*!
A function object to construct a `Segment_2`.

Provides:

`Segment_2 operator()(Point_2 p,Point_2 q)`,

which constructs a segment from two points.
*/
typedef unspecified_type Construct_segment_2;

/*!
A function object to construct a `Triangle_2`.

Provides:

`Triangle_2 operator()(Point_2 p,Point_2 q,Point_2 r )`,

which constructs a triangle from three points.
*/
typedef unspecified_type Construct_triangle_2;

/*!
A function object to compare the x-coordinate of two points.

Provides the operator:

`bool operator()(Point p, Point q)`

which returns `true` if `p` is before `q`
according to the \f$ x\f$-ordering of points.
*/
typedef unspecified_type Less_x_2;

/*!
A function object to compare the y-coordinate of two points.

Provides the operator:

`bool operator()(Point p, Point q)`

which returns `true` if `p` is before `q`
according to the \f$ y\f$-ordering of points.
*/
typedef unspecified_type Less_y_2;

/*!
A function object to compare the x-coordinate of two points.

Provides the operator:

`Comparison_result operator()(Point p, Point q)`

which returns
`SMALLER, EQUAL` or `LARGER`
according to the
\f$ x\f$-ordering of points `p` and `q`.
*/
typedef unspecified_type Compare_x_2;

/*!
A function object to compare the y-coordinate of two points.
Provides the operator:

`Comparison_result operator()(Point p, Point q)`

which returns
(`SMALLER, EQUAL` or `LARGER`)
according to the
\f$ y\f$-ordering of points `p` and `q`.
*/
typedef unspecified_type Compare_y_2;

/*!
A function object to compute the orientation of three points.

Provides the operator:

`Orientation operator()(Point p, Point q, Point r)`

which returns `LEFT_TURN`, `RIGHT_TURN` or `COLLINEAR`
depending on \f$ r\f$ being, with respect to
the oriented line `pq`,
on the left side , on the right side or on the line.
*/
typedef unspecified_type Orientation_2;

/*!
A function object to perform the incircle test for four points.

Provides the operator:

`Oriented_side operator()(Point p, Point q, Point r, Point s)`
which takes four points \f$ p, q, r, s\f$ as arguments and returns
`ON_POSITIVE_SIDE`, `ON_NEGATIVE_SIDE` or,
`ON_ORIENTED_BOUNDARY` according to the position of points `s`
with respect to the oriented circle through through \f$ p,q\f$
and \f$ r\f$.
This type is required only if the function
`side_of_oriented_circle(Face_handle f, Point p)` is
called.
*/
typedef unspecified_type Side_of_oriented_circle_2;

/*!
A function object to compute the circumcentr of three points.
Provides the operator:

`Point operator()(Point p, Point q, Point r)`

which returns
the circumcenter of the three points `p, q` and `r`.
This type is required only if the function
`Point circumcenter(Face_handle f)`is called.
*/
typedef unspecified_type Construct_circumcenter_2;

/// @}

/// \name Creation
/// Only a default constructor, copy constructor and an assignment
/// operator are required. Note that further constructors can be
/// provided.
/// @{

/*!
default constructor.
*/
TriangulationTraits_2();

/*!
Copy constructor
*/
TriangulationTraits_2(TriangulationTraits_2 gtr);

/*!
Assignment operator.
*/
TriangulationTraits_2 operator=(TriangulationTraits_2 gtr);

/// @}

/// \name Predicate Functions
/// The following functions give access to the predicate and
/// constructor objects.
/// @{

/*!

*/
Construct_point_2 construct_point_2_object();

/*!

*/
Construct_segment_2 construct_segment_2_object();

/*!

*/
Construct_triangle_2 construct_triangle_2_object();

/*!

*/
Compare_x_2 compare_x_2_object();

/*!

*/
Compare_y_2 compare_y_2_object();

/*!

*/
Orientation_2 orientation_2_object();

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

/// @}

}; /* end TriangulationTraits_2 */

