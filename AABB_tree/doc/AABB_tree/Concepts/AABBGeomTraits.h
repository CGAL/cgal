
/*!
\ingroup PkgAABBTreeConcepts
\cgalConcept

The concept `AABBGeomTraits` defines the requirements for the first template parameter of the class
`CGAL::AABB_traits<AABBGeomTraits, AABBPrimitive>`. It provides predicates and constructors to detect
and compute intersections between query objects and the primitives stored in the AABB tree.
In addition, it contains predicates and constructors to compute distances between a point query
and the primitives stored in the AABB tree.

\cgalRefines{SearchGeomTraits_3}

\cgalHasModelsBegin
\cgalHasModelsBare{All models of the concept `Kernel`}
\cgalHasModelsEnd

\sa `CGAL::AABB_traits<AABBGeomTraits,AABBPrimitive>`
\sa `CGAL::AABB_tree<AABBTraits>`
\sa `AABBPrimitive`

*/

class AABBGeomTraits {
public:

/// \name Types
/// @{

/*!
A functor object to detect intersections between two geometric objects.
Provides the following operators:

`bool operator()(Query, Bbox_3)`,

`bool operator()(Query, Primitive::Datum)`,

`bool operator()(Sphere_3, Bbox_3)`.

The operator returns `true` iff there exists a non-empty intersection.
*/
typedef unspecified_type Do_intersect_3;

/*!
A functor object to construct the intersection between two geometric objects.

Provides the operator:

`return_type operator()(const Query& q, const Primitive::Datum& d)`,

which computes the intersection between `q` and `d`. The type of the returned object
must be a `std::optional` of a `std::variant` of the possible intersection types.
*/
typedef unspecified_type Intersect_3;

/*!
A functor object to construct the sphere centered at one point and passing through another one.
Provides the operator:

`Sphere_3 operator()(const Point_3& p, const FT & sr)`,

which returns the sphere centered at `p` with `sr` as squared radius.
*/
typedef unspecified_type Construct_sphere_3;

/*!
A functor object to compute the point on a geometric primitive which is closest from a query.
Provides the operator:

`Point_3 operator()(const Primitive::Datum& d, const Point_3& p)`,

which returns the point on `d` that is closest to `p`.
*/
typedef unspecified_type Construct_projected_point_3;

/*!
A functor object to compare the distance of two points wrt a third one. Provides the operator:

`CGAL::Comparison_result operator()(const Point_3& p1, const Point_3& p2, const Point_3& p3)`,

which compares the distance between `p1` and `p2`, and between `p2` and `p3`.
*/
typedef unspecified_type Compare_distance_3;

/*!
A functor object to compute the squared radius of a sphere.
Provides the operator:

`FT operator()(const Sphere_3& s),`

which returns the squared radius of `s`.
*/
typedef unspecified_type Compute_squared_radius_3;

/*!
A functor object to compute the squared distance between two points. Provides the operator:

`FT operator()(const Point_3& p, const Point_3& q),`

which returns the squared distance between `p` and `q`.
*/
typedef unspecified_type Compute_squared_distance_3;

/*!
A functor object to compare the x-coordinates of two points. Provides the operator:

`bool operator()(const Point_3& p, const Point_3& q)`,

 which returns `true` iff the x-coordinate of `p` is smaller than the x-coordinate of `q`.
*/
typedef unspecified_type Less_x_3;

/*!
A functor object to compare the y-coordinates of two points. Provides the operator:

`bool operator()(const Point_3& p, const Point_3& q)`,

which returns `true` iff the y-coordinate of `p` is smaller than the y-coordinate of `q`.
*/
typedef unspecified_type Less_y_3;

/*!
A functor object to compare the z-coordinates of two points. Provides the operator:

`bool operator()(const Point_3& p, const Point_3& q)`,

which returns `true` iff the z-coordinate of `p` is smaller than the z-coordinate of `q`.
*/
typedef unspecified_type Less_z_3;

/*!
A functor object to compare two points. Provides the operator:

`bool operator()(const Point_3& p, const Point_3& q)`,

which returns `true` iff `p` is equal to `q`.
*/
typedef unspecified_type Equal_3;

/// @}

/// \name Operations
/// @{

/*!
returns the intersection detection predicate.
*/
Do_intersect_3 do_intersect_3_object();

/*!
returns the intersection constructor.
*/
Intersect_3 intersect_3_object();

/*!
returns the sphere constructor.
*/
Construct_sphere_3 construct_sphere_3_object();

/*!
returns the closest point constructor.
*/
Construct_projected_point_3 construct_projected_point_3_object();

/*!
returns the compare distance predicate.
*/
Compare_distance_3 compare_distance_3_object();

/*!
returns the squared radius functor.
*/
Compute_squared_radius_3 compute_squared_radius_3_object();

/*!
returns the squared distance functor.
*/
Compute_squared_distance_3 compute_squared_distance_3_object();

/*!
returns the `Less_x_3` predicate.
*/
Less_x_3 less_x_3_object();

/*!
returns the `Less_y_3` predicate.
*/
Less_y_3 less_y_3_object();

/*!
returns the `Less_z_3` predicate.
*/
Less_z_3 less_z_3_object();

/*!
returns the equal predicate.
*/
Equal_3 equal_3_object();

/// @}

}; /* end AABBGeomTraits */

