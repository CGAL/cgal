/*!
\ingroup PkgConvexHull3Concepts
\cgalConcept

Requirements of the traits class of the function `CGAL::convex_hull_3()`. 

\cgalHasModel `CGAL::Convex_hull_traits_3` 
\cgalHasModel All models of `Kernel`

*/
class ConvexHullTraits_3 {
public:

/// \name Types 
/// @{

/*!
The point type on which the convex hull algorithm operates 
*/ 
typedef unspecified_type Point_3; 

/*!
a 3D plane 
*/ 
typedef unspecified_type Plane_3; 

/*!
a 3D segment 
*/ 
typedef unspecified_type Segment_3; 

/*!
a 3D triangle 
*/ 
typedef unspecified_type Triangle_3; 

/*!
Function object type that provides 
`Plane_3 operator()(Point_3 p, Point_3 q, Point_3 r)`, which constructs 
and returns a plane passing through `p`, `q`, and `r` and oriented 
in a positive sense when seen from the positive side of the plane. 
*/ 
typedef unspecified_type Construct_plane_3; 

/*!
Function object type that provides 
`Segment_3 operator()(Point_3 p, Point_3 q)`, which constructs and 
returns the segment with source `p` and target `q`. 
*/ 
typedef unspecified_type Construct_segment_3; 

/*!
Function object type that provides 
`Triangle_3 operator()(Point_3 p, Point_3 q, Point_3 r)`, which 
constructs and returns the triangle with vertices `p`, `q`, and 
`r`. 
*/ 
typedef unspecified_type Construct_triangle_3;  

/*!
Predicate object type that provides 
`bool operator()(Point_3 p, Point_3 q)`, which determines 
if points `p` and `q` are equal. 
*/ 
typedef unspecified_type Equal_3; 

/*!
Predicate object type that provides 
`bool operator()(Point_3 p, Point_3 q, Point_3 r)`, which determines 
if points `p`, `q` and `r` are collinear. 
*/ 
typedef unspecified_type Collinear_3; 

/*!
Predicate object type that provides 
`bool operator()(Point_3 p, Point_3 q, Point_3 r, Point_3 s)`, which 
determines if points `p`, `q`, `r`, and `s` are coplanar. 
*/ 
typedef unspecified_type Coplanar_3; 

/*!
Predicate object type that provides 
`bool operator()(Plane_3 h, Point_3 q)`, which determines if the point 
`q` is on the positive side of the halfspace `h`. 
*/ 
typedef unspecified_type Has_on_positive_side_3; 

/*!
Predicate object type that provides 
a constructor taking a single `Point_3` object and 
`bool operator()(Point_3 q, Point_3 r)`, which returns true iff the 
distance from `q` to `p` is smaller than the distance from 
`r` to `p`, where `p` is the point passed to the object 
at construction. 
*/ 
typedef unspecified_type Less_distance_to_point_3; 

/*!
Predicate object type that 
provides `bool operator()(Plane_3 p, Point_3 q, Point_3 r)`, which 
returns true iff the signed distance from `q` to `p` is smaller 
than the signed distance from `r` to `p` 
*/ 
typedef unspecified_type Less_signed_distance_to_plane_3; 

/*!
A traits class providing the requirements of the template parameter `Traits` of
the 2D convex hull function `CGAL::ch_bykat()` such that `Traits::Point_2`
is `Point_3`, and the 2D points considered in the algorithm are the projections
of the 3D points in the `xy`-plane.
If this type is not available, the function `CGAL::convex_hull_3()` will
automatically use `CGAL::Projection_traits_xy< CGAL::Kernel_traits<Point_3>::%Kernel >.`
*/
typedef unspecified_type Traits_xy_3;

/*!
Same as above but in the `yz`-plane
*/
typedef unspecified_type Traits_yz_3;

/*!
Same as above but in the `xz`-plane
*/
typedef unspecified_type Traits_xz_3;

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
