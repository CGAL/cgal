
/*!
\ingroup PkgAABBTreeConcepts
\cgalConcept

The concept `AABBGeomTraits_2` defines the requirements for the first template parameter of the class
`CGAL::AABB_traits_2<AABBGeomTraits_2, AABBPrimitive>`. It provides predicates and constructors to detect
and compute intersections between query objects and the primitives stored in the AABB tree.
In addition, it contains predicates and constructors to compute distances between a point query
and the primitives stored in the AABB tree.

\cgalRefines{SearchGeomTraits_2}

\cgalHasModelsBegin
\cgalHasModelsBare{All models of the concept `Kernel`}
\cgalHasModelsEnd

\sa `CGAL::AABB_traits_2<AABBGeomTraits_2,AABBPrimitive>`
\sa `CGAL::AABB_tree<AABBTraits>`
\sa `AABBPrimitive`

*/

class AABBGeomTraits_2 {
public:

/// \name Types
/// @{

/*!
A functor object to detect intersections between two geometric objects.
Provides the following operators:

`bool operator()(const Query& q, const Bbox_2& b)`,

`bool operator()(const Query& q, const Primitive::Datum& d)`,

`bool operator()(const Circle_2& c, const Bbox_2& b)`.

The operator returns `true` iff there is an intersection.
*/
typedef unspecified_type Do_intersect_2;

/*!
A functor object to construct the intersection between two geometric objects.

Provides the operator:

`return_type operator()(const Query& q, const Primitive::Datum& d)`,

which computes the intersection between `q` and `d`. The type of the returned object
must be a `std::optional` of a `std::variant` of the possible intersection types.
*/
typedef unspecified_type Intersect_2;

/*!
A functor object to construct the circle specified by its center and squared radius.
Provides the operator:

`Circle_2 operator()(const Point_2& p, const FT & sr)`,

which returns the circle centered at `p` with `sr` as squared radius.
*/
typedef unspecified_type Construct_circle_2;

/*!
A functor object to compute the point on a geometric primitive which is closest from a query point.
Provides the operator:

`Point_2 operator()(const Primitive::Datum& d, const Point_2& p)`,

which returns the point on `d` that is closest to `p`.
*/
typedef unspecified_type Construct_projected_point_2;

/*!
A functor object to compare the distance of two points wrt a third one. Provides the operator:

`CGAL::Comparison_result operator()(const Point_2& p1, const Point_2& p2, const Point_2& p3)`,

which compares the distance between `p1` and `p2`, to the distance between `p1` and `p3`.
*/
typedef unspecified_type Compare_distance_2;


/*!
A functor object to compute the squared distance between two points. Provides the operator:

`FT operator()(const Point_2& p, const Point_2& q),`

which returns the squared distance between `p` and `q`.
*/
typedef unspecified_type Compute_squared_distance_2;

/*!
A functor object to compare the x-coordinates of two points. Provides the operator:

`bool operator()(const Point_2& p, const Point_2& q)`,

 which returns `true` iff the x-coordinate of `p` is smaller than the x-coordinate of `q`.
*/
typedef unspecified_type Less_x_2;

/*!
A functor object to compare the y-coordinates of two points. Provides the operator:

`bool operator()(const Point_2& p, const Point_2& q)`,

which returns `true` iff the y-coordinate of `p` is smaller than the y-coordinate of `q`.
*/
typedef unspecified_type Less_y_2;


/*!
A functor object to compare two points. Provides the operator:

`bool operator()(const Point_2& p, const Point_2& q)`,

which returns `true` iff `p` is equal to `q`.
*/
typedef unspecified_type Equal_2;

/// @}

/// \name Operations
/// @{

/*!
returns the intersection detection predicate.
*/
Do_intersect_2 do_intersect_2_object();

/*!
returns the intersection constructor.
*/
Intersect_2 intersect_2_object();

/*!
returns the circle constructor.
*/
Construct_circle_2 construct_circle_2_object();

/*!
returns the closest point constructor.
*/
Construct_projected_point_2 construct_projected_point_2_object();

/*!
returns the compare distance predicate.
*/
Compare_distance_2 compare_distance_2_object();


/*!
returns the squared distance functor.
*/
Compute_squared_distance_2 compute_squared_distance_2_object();

/*!
returns the `Less_x_2` predicate.
*/
Less_x_2 less_x_2_object();

/*!
returns the `Less_y_2` predicate.
*/
Less_y_2 less_y_2_object();


/*!
returns the equal predicate.
*/
Equal_2 equal_2_object();

/// @}

}; /* end AABBGeomTraits_2 */
