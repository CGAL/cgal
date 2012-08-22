
/*!
\ingroup PkgTriangulation2Concepts
\cgalconcept

The concept `ConstrainedTriangulationTraits_2` defines the requirements for the geometric 
traits class of a constrained triangulation 
( `Constrained_Triangulation_2<Traits,Tds,Itag>`) 
that supports intersections of input constraints (i. e. 
when the template parameter `Itag` is instantiated 
by one of the tag classes `Exact_intersections_tag` or 
`Exact_predicates_tag`). This concept refines the concept 
`TriangulationTraits_2`, adding requirements for function objects 
to compute the intersection points of two constraints. 
When `Exact_predicates_tag` is used, the 
traits class is 
also required to provide additional types 
to compute the squared distance between a point and a line 

\refines ::TriangulationTraits_2 

\hasModel All \cgal Kernels
\hasModel CGAL::Projection_traits_xy_3<K> 
\hasModel CGAL::Projection_traits_yz_3<K> 
\hasModel CGAL::Projection_traits_zx_3<K> 

\sa `TriangulationTraits_2` 
\sa `ConstrainedDelaunayTriangulationTraits_2` 
\sa `CGAL:Constrained_Triangulation_2<Traits,Tds,Itag>` 

*/

class ConstrainedTriangulationTraits_2 {
public:

/// \name Types 
/// @{

/*! 
A function object whose operator() computes the intersection of two segments: 

`Object_2 operator()(Segment_2 s1, Segment_2 s2);` 
Returns the intersection of `s1` and `s2`. 
*/ 
typedef Hidden_type Intersect_2; 

/// \name Types requires with ::Exact_predicates_tag
/// When the constrained triangulation is instantiated with the intersection tag `Exact_predicates_tag`, the used algorithm needs to be able to compare some distances between points and lines and the following types are further required.
/// @{

/*! 
A number type supporting the comparison operator 
`<`. 
*/ 
typedef Hidden_type RT; 

/*! 
The line type. 
*/ 
typedef Hidden_type Line_2; 

/*! 
A function object whose operator() 
constructs a line from two points : 

`Line_2 operator()(Point_2 p1, Point_2 p2)`. 
*/ 
typedef Hidden_type Construct_line_2; 

/*! 
A function object with an 
operator() designed to compute the squared distance between 
a line and a point : 

`RT operator()(Line_2 l, Point_2 p);` Return the squared distance 
between `p` and `l`. 
*/ 
typedef Hidden_type Compute_squared_distance_2; 

/// @} 

/// \name Access to constructor object 
/// @{

/*! 

*/ 
Intersect_2 intersect_2_object(); 

/*! 
required when 
the intersection tag is `Exact_predicates_tag`. 
*/ 
Construct_line_2 construct_line_2_object(); 

/*! 
required when 
the intersection tag is `Exact_predicates_tag`. 
*/ 
Compute_squared_distance_2 
compute_squared_distance_2_object(); 

/// @}

}; /* end ConstrainedTriangulationTraits_2 */

