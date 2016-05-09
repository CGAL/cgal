
/*!
\ingroup PkgAABB_treeConcepts
\cgalConcept

The concept `AABBGeomTraits` defines the requirements for the first template parameter of the class `CGAL::AABB_traits<AABBGeomTraits, AABBPrimitive>`. It provides predicates and constructors to detect and compute intersections between query objects and the primitives stored in the AABB tree. In addition, it contains predicates and constructors to compute distances between a point query and the primitives stored in the AABB tree. 

\cgalHasModel Any 3D Kernel is a model of this traits concept.

\sa `CGAL::AABB_traits<AABBGeomTraits,AABBPrimitive>`

*/

class AABBGeomTraits {
public:

/// \name Types 
/// @{

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
This functor must support the result_of protocol, that is the return 
type of the `operator()(A, B)` is `CGAL::cpp11::result<Intersect_3(A,B)>`.

Provides the operators: 
`CGAL::cpp11::result<Intersect_3(A,B)> operator()(const A& a, const B& b);` 
where `A` and `B` are any relevant types among `Ray_3`, `Segment_3`, `Line_3`, 
`Triangle_3`, `Plane_3` and `Bbox_3`. 
Relevant herein means that a line primitive (ray, segment, line) is tested 
against a planar or solid primitive (plane, triangle, box). 
A model of `Kernel::Intersect_3` fulfills those requirements. 
*/ 
typedef unspecified_type Intersect_3; 

/*!
A functor object to construct the sphere centered at one point and passing through another one. Provides the operator: 
`Sphere_3 operator()(const Point_3& p, const Point_3 & q);` which returns the sphere centered at `p` and passing through `q`. 
*/ 
typedef unspecified_type Construct_sphere_3; 

/*!
A functor object to compute the point on a geometric primitive which is closest from a query. Provides the operator: 
`Point_3 operator()(const Point_3& p, const Type_2& type_2);` where `Type_2` is any type among `Segment_3` and `Triangle_3`. The operator returns the point on `type_2` which is closest to `p`. 
*/ 
typedef unspecified_type Compute_closest_point_3; 

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
`FT operator()(const Point_3& p, const Point_3& q);}` which returns the squared distance between \a p and \a q.
*/
typedef unspecified_type Compute_squared_distance_3;


/// @} 

/// \name Operations 
/// @{

/*!
Returns the intersection detection functor. 
*/ 
Do_intersect_3 do_intersect_3_object(); 

/*!
Returns the intersection constructor. 
*/ 
Intersect_3 intersect_3_object(); 

/*!
Returns the sphere constructor.
*/ 
Construct_sphere_3 construct_sphere_3_object(); 

/*!
Returns the closest point constructor. 
*/ 
Compute_closest_point_3 compute_closest_point_3_object(); 

/*!
Returns the closest point constructor. 
*/ 
Has_on_bounded_side_3 has_on_bounded_side_3_object(); 

/*!
Returns the squared radius functor. 
*/ 
Compute_squared_radius_3 compute_squared_radius_3_object(); 

/*!
Returns the squared distance functor.
*/
Compute_squared_distance_3 compute_squared_distance_3_object();

/// @}

}; /* end AABBGeomTraits */

