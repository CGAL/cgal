/*!
\ingroup PkgConvexHull3Concepts
\cgalconcept

Requirements of the traits class to be used with the function 
`convex_hull_3`. 

\hasModel CGAL::Convex_hull_traits_3<R> 
\hasModel All kernels of CGAL 

*/
class ConvexHullTraits_3 {
public:

/// \name Types 
/// @{

/*! 
The point type on which the convex hull algorithm operates 
*/ 
typedef Hidden_type Point_3; 

/*! 
a 3D plane 
*/ 
typedef Hidden_type Plane_3; 

/*! 
a 3D segment 
*/ 
typedef Hidden_type Segment_3; 

/*! 
a 3D triangle 
*/ 
typedef Hidden_type Triangle_3; 

/*! 
a 3D vector 
*/ 
typedef Hidden_type Vector_3; 

/*! 
Function object type that provides 
`Plane_3 operator()(Point_3 p, Point_3 q, Point_3 r)`, which constructs 
and returns a plane passing through `p`, `q`, and `r` and oriented 
in a positive sense when seen from the positive side of the plane. 
*/ 
typedef Hidden_type Construct_plane_3; 

/*! 
Function object type that provides 
`Segment_3 operator()(Point_3 p, Point_3 q)`, which constructs and 
returns the segment with source `p` and target `q`. 
*/ 
typedef Hidden_type Construct_segment_3; 

/*! 
Function object type that provides 
`Triangle_3 operator()(Point_3 p, Point_3 q, Point_3 r)`, which 
constructs and returns the triangle with vertices `p`, `q`, and 
`r` 
*/ 
typedef Hidden_type Construct_triangle_3; 

/*! 
Function object type that provides 
`Vector_3 operator()(Point_3 p, Point_3 q)`, which constructs and 
returns the vector `q`-`p` 
*/ 
typedef Hidden_type Construct_vector_3; 

/*! 
Predicate object type that provides 
`bool operator()(Point_3 p, Point_3 q)`, which determines 
if points `p` and `q` are equal or not 
*/ 
typedef Hidden_type Equal_3; 

/*! 
Predicate object type that provides 
`bool operator()(Point_3 p, Point_3 q, Point_3 r)`, which determines 
if points `p`, `q` and `r` are collinear or not 
*/ 
typedef Hidden_type Collinear_3; 

/*! 
Predicate object type that provides 
`bool operator()(Point_3 p, Point_3 q, Point_3 r, Point_3 s)`, which 
determines if points `p`, `q`, `r`, and `s` are coplanar 
or not 
*/ 
typedef Hidden_type Coplanar_3; 

/*! 
Predicate object type that provides 
`bool operator()(Plane_3 h, Point_3 q)`, which determines of the point 
`q` is on the positive side of the halfspace `h` 
*/ 
typedef Hidden_type Has_on_positive_side_3; 

/*! 
Predicate object type that provides 
a constructor taking a single `Point_3` object and 
`bool operator()(Point_3 q, Point_3 r)`, which returns true iff the 
distance from `q` to `p` is smaller than the distance from 
`r` to `p`, where `p` is the point passed to the object 
at construction. 
*/ 
typedef Hidden_type Less_distance_to_point_3; 

/*! 
Predicate object type that 
provides `bool operator()(Plane_3 p, Point_3 q, Point_3 r)`, which 
returns true iff the signed distance from `q` to `p` is smaller 
than the signed distance from `r` to `p` 
*/ 
typedef Hidden_type Less_signed_distance_to_plane_3; 

/// @} 

/// \name Creation 
/// Only a copy constructor is required. 
/// @{

/*! 

*/ 
ConvexHullTraits_3(ConvexHullTraits_3& ch); 

/// @} 

/*! \name Operations 
For each of the above function and predicate object types, 
`Func_obj_type`, a function must exist with the name 
`func_obj_type_object` that creates an instance of the function or 
predicate object type. For example: 
*/
/// @{

/*! 

*/ 
Construct_plane_3 construct_plane_3_object(); 

/// @}

}; /* end ConvexHullTraits_3 */
