


/*!
\ingroup PkgMesh2Concepts
\cgalconcept


The concept `DelaunayMeshTraits_2` refines the concept 
`ConformingDelaunayTriangulationTraits_2`. It provides a construction 
object `Construct_circumcenter_2`. 

\refines ::ConformingDelaunayTriangulationTraits_2 

\hasModel Any model of the `Kernel` concept. In particular, all \cgal kernels. 


*/

class DelaunayMeshTraits_2 {
public:


/// \name Types 
/// @{

/*! 
Constructor object. Must provide 
an operator `Point_2 operator()(Point_2 p, Point_2 q, Point_2 r);` 
that constructs the center of the circle passing through the points `p`, 
`q`, and `r`. 
\pre `p`, `q`, and `r` are not collinear. 
*/ 
typedef Hidden_type Construct_circumcenter_2; 


/*! 
Constructor object. Must provide an operator 
`FT operator()(Point_2 p, Point_2 q, Point_2 r);` 
that computes the signed area of the triangle defined by 
the points `p`, `q`, and `r`. 
*/ 
typedef Hidden_type Compute_area_2; 

/// @} 


/// \name Access to predicate and constructor objects 
/// @{

/*! 

*/ 
Construct_circumcenter_2 construct_circumcenter_2_object(); 





/*! 

*/ 
Compute_area_2 compute_area_2_object(); 





/// @}

}; /* end DelaunayMeshTraits_2 */

