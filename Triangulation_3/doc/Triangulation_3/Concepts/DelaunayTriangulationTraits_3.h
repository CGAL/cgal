
/*!
\ingroup PkgTriangulation3Concepts
\cgalConcept

The concept `DelaunayTriangulationTraits_3` is the first template parameter of the class 
`Delaunay_triangulation_3`. It defines the geometric objects (points, 
segments...) forming the triangulation together with a few geometric 
predicates and constructions on these objects. 

\cgalRefines `TriangulationTraits_3` 

\cgalHasModel `CGAL::Exact_predicates_inexact_constructions_kernel` (recommended) 
\cgalHasModel `CGAL::Exact_predicates_exact_constructions_kernel` (recommended for Voronoi) 
\cgalHasModel CGAL::Filtered_kernel 
\cgalHasModel CGAL::Cartesian 
\cgalHasModel CGAL::Simple_cartesian 
\cgalHasModel CGAL::Homogeneous 
\cgalHasModel CGAL::Simple_homogeneous 

In addition to the requirements described for the traits class of 
`CGAL::Triangulation_3`, the geometric traits class of a 
  Delaunay triangulation must fulfill the following requirements: 

*/

class DelaunayTriangulationTraits_3 {
public:

/// \name Types 

/// @{

/*!
The line type. 
*/ 
typedef unspecified_type Line_3; 

/*!
The object type. 
*/ 
typedef unspecified_type Object_3; 

/*!
The ray type. 
*/ 
typedef unspecified_type Ray_3; 

/*!
A predicate object that must provide the function operator 

`Bounded_side operator()(Point p, Point q, Point r, Point s)`, 

which determines the bounded side of the circle defined 
by `p`, `q`, and `r` on which `s` lies. 
\pre `p`, `q`, `r`, and `s` are coplanar and `p`, `q`, and `r` are not collinear. 
*/ 
typedef unspecified_type Coplanar_side_of_bounded_circle_3; 

/*!
A predicate object that must provide the function operator 

`Oriented_side operator()(Point p, Point q, Point r, Point s, Point t)`, 

which determines on which side of the oriented sphere circumscribing 
`p, q, r, s` the point `t` lies. 
*/ 
typedef unspecified_type Side_of_oriented_sphere_3; 

/*!
A predicate object that must provide the function operator 

`Comparison_result operator()(Point p, Point q, Point r)`, 

which compares the distance between `p` and `q` to the distance 
between `p` and `r`. 

It is only needed when using the `Fast_location` policy or the 
`nearest_vertex` function. 
*/ 
typedef unspecified_type Compare_distance_3; 


/// @}

/*! \name
 In addition, only when the dual operations are used, the traits class must provide the following constructor objects:
*/
/// @{

/*!
A constructor object that must provide the function operator 

`Point_3 operator()(Point_3 p, Point_3 q, Point_3 r, Point_3 s)`, 

which constructs the circumcenter of four points. 
\pre `p`, `q`, `r` and `s` must be non coplanar. 

It must also provide the function operator 

`Point_3 operator()(Point_3 p, Point_3 q, Point_3 r)`, 

which constructs the circumcenter of three points. 
\pre `p`, `q` and `r` must be non collinear. 
*/ 
typedef unspecified_type Construct_circumcenter_3; 

/*!
A constructor object that must provide the function operators 

`Object_3 operator()(Point_3 p)`, 

`Object_3 operator()(Segment_3 s)` and 

`Object_3 operator()(Ray_3 r)` 

that construct an object respectively from a point, a segment and a ray. 
*/ 
typedef unspecified_type Construct_object_3; 

/*!
A constructor object that must provide the function operator 

`Line_3 operator()(Point_3 p1, Point_3 p2, Point_3 p3)`, 

which constructs the line which is at the same distance from the three points. 
\pre `p1`, `p2` and `p3` must be non collinear. 
*/ 
typedef unspecified_type Construct_equidistant_line_3; 

/*!
A constructor object that must provide the function operator 

`Ray_3 operator()(Point_3 p, Line_3 l)`, 

which constructs the ray starting at `p` with direction given by `l`. 
*/ 
typedef unspecified_type Construct_ray_3; 

/// @} 

/// \name Operations 
/// The following functions give access to the predicate and construction objects:
/// @{

/*!

*/ 
Coplanar_side_of_bounded_circle_3 coplanar_side_of_bounded_circle_3_object(); 

/*!

*/ 
Side_of_oriented_sphere_3 side_of_oriented_sphere_3_object(); 


/// @}

/*! \name
When using the `Fast_location` policy or the `CGAL::Delaunay_triangulation_3::nearest_vertex()` function, the traits must provide:
*/
/// @{

/*!

*/ 
Compare_distance_3 compare_distance_3_object(); 



/// @}

/*! \name
The following functions must be provided only if the methods of `Delaunay_triangulation_3` returning elements of the Voronoi diagram are instantiated:
*/
/// @{

/*!

*/ 
Construct_circumcenter_3 construct_circumcenter_3_object(); 

/*!

*/ 
Construct_object_3 construct_object_3_object(); 

/*!

*/ 
Construct_equidistant_line_3 construct_equidistant_line_object(); 

/*!

*/ 
Construct_ray_3 construct_ray_3_object(); 

/// @}

}; /* end DelaunayTriangulationTraits_3 */

