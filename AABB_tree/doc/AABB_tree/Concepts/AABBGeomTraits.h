
/*!
\ingroup PkgAABBTreeConcepts
\cgalConcept

The concept `AABBGeomTraits` defines the requirements for the first template parameter of the class `CGAL::AABB_traits<AABBGeomTraits, AABBPrimitive>`. It provides predicates and constructors to detect and compute intersections between query objects and the primitives stored in the AABB tree. In addition, it contains predicates and constructors to compute distances between a point query and the primitives stored in the AABB tree.

\cgalRefines `SearchGeomTraits_3`

\cgalHasModel Any 3D Kernel is a model of this traits concept.

\sa `CGAL::AABB_traits<AABBGeomTraits,AABBPrimitive>`

*/

class AABBGeomTraits {
public:

/// \name Types
/// @{

/*!
A number type model of `Field`.
*/
typedef unspecified_type FT;

/*!
Sphere type, that should be consistent with the distance function chosen for the distance queries, namely the `Squared_distance_3` functor.
*/
typedef unspecified_type Sphere_3;

/*!
Point type.
*/
typedef unspecified_type Point_3;

/*!
A functor object to detect intersections between two geometric objects.
Provides the operators:
`bool operator()(const Type_1& type_1, const Type_2& type_2);`
where `Type_1` and `Type_2` are relevant types
among `Ray_3`, `Segment_3`, `Line_3`, `Triangle_3`, `Plane_3` and `Bbox_3`. Relevant herein means that a line primitive (ray, segment, line) is tested against a planar or solid primitive (plane, triangle, box), and a solid primitive is tested against another solid primitive (box against box). The operator returns `true` iff `type_1` and `type_2` have a non empty intersection.
*/
typedef unspecified_type Do_intersect_3;

/*!
A functor object to construct the intersection between two geometric objects.

Provides the operators:
`decltype(auto) operator()(const A& a, const B& b);`
where `A` and `B` are any relevant types among `Ray_3`, `Segment_3`, `Line_3`,
`Triangle_3`, `Plane_3` and `Bbox_3`.
Relevant herein means that a line primitive (ray, segment, line) is tested
against a planar or solid primitive (plane, triangle, box).
A model of `Kernel::Intersect_3` fulfills those requirements.
*/
typedef unspecified_type Intersect_3;

/*!
A functor object to construct the sphere centered at one point and passing through another one. Provides the operator:
- `Sphere_3 operator()(const Point_3& p, const FT & sr)` which returns the sphere centered at `p` with `sr` as squared radius.
*/
typedef unspecified_type Construct_sphere_3;

/*!
A functor object to compute the point on a geometric primitive which is closest from a query. Provides the operator:
`Point_3 operator()(const Type_2& type_2, const Point_3& p);` where `Type_2` can be any of the following types : `Segment_3`, `Ray_3`, or `Triangle_3`.
The operator returns the point on `type_2` which is closest to `p`.
*/
typedef unspecified_type Construct_projected_point_3;

/*!
A functor object to compare the distance of two points wrt a third one.
Provides the operator:
`CGAL::Comparision_result operator()(const Point_3& p1, const Point_3& p2, const Point_3& p3)`. The operator compare the distance between `p1 and `p2`, and between `p2` and `p3`.
*/
typedef unspecified_type Compare_distance_3;

/*!
A functor object to detect if a point lies inside a sphere or not.
Provides the operator:
`bool operator()(const Sphere_3& s, const Point_3& p);` which returns `true` iff the closed volume bounded by `s` contains `p`.
*/
typedef unspecified_type Has_on_bounded_side_3;

/*!
A functor object to compute the squared radius of a sphere. Provides the operator:
`FT operator()(const Sphere_3& s);` which returns the squared radius of `s`.
*/
typedef unspecified_type Compute_squared_radius_3;

/*!
A functor object to compute the squared distance between two points. Provides the operator:
`FT operator()(const Point_3& p, const Point_3& q);}` which returns the squared distance between `p` and `q`.
*/
typedef unspecified_type Compute_squared_distance_3;

/*!
A functor object to compare the x-coordinates of two points. Provides the operator:
`bool operator()(const Point_3& p, const Point_3& q);}` which returns `true` if the x-coordinate of `p` is smaller
than the x-coordinate of `q`.
*/
typedef unspecified_type Less_x_3;

/*!
A functor object to compare the y-coordinates of two points. Provides the operator:
`bool operator()(const Point_3& p, const Point_3& q);}` which returns `true` if the y-coordinate of `p` is smaller
than the y-coordinate of `q`.
*/
typedef unspecified_type Less_y_3;

/*!
A functor object to compare the z-coordinates of two points. Provides the operator:
`bool operator()(const Point_3& p, const Point_3& q);}` which returns `true` if the z-coordinate of `p` is smaller
than the z-coordinate of `q`.
*/
typedef unspecified_type Less_z_3;

/*!
A functor object to compare two points. Provides the operator:
`bool operator()(const Point_3& p, const Point_3& q);}` which returns `true` if `p` is equal to `q`.
*/
typedef unspecified_type Equal_3;



/// @}

/// \name Operations
/// @{

/*!
returns the intersection detection functor.
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
returns the compare distance constructor.
*/
Compare_distance_3 compare_distance_3_object();

/*!
returns the closest point constructor.
*/
Has_on_bounded_side_3 has_on_bounded_side_3_object();

/*!
returns the squared radius functor.
*/
Compute_squared_radius_3 compute_squared_radius_3_object();

/*!
returns the squared distance functor.
*/
Compute_squared_distance_3 compute_squared_distance_3_object();

/*!
returns the `Less_x_3` functor.
*/
Less_x_3 less_x_3_object();

/*!
returns the `Less_y_3` functor.
*/
Less_y_3 less_y_3_object();

/*!
returns the `Less_z_3` functor.
*/
Less_z_3 less_z_3_object();

/*!
returns the equal functor.
*/
Equal_3 equal_3_object();

/// @}

}; /* end AABBGeomTraits */

