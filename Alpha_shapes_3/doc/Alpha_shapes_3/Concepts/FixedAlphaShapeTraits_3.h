
/*!
\ingroup PkgAlphaShapes3Concepts
\cgalconcept

The concept `FixedAlphaShapeTraits_3` describes the requirements 
for the geometric traits class 
of the underlying Delaunay triangulation of a basic alpha shape with a fixed value alpha. 

\refines ::DelaunayTriangulationTraits_3 

In addition to the requirements described in the concept
`DelaunayTriangulationTraits_3`, the geometric traits class of a
Delaunay triangulation plugged in a basic alpha shape with fixed alpha
value provides the following.

\hasModel All CGAL kernels. 
\hasModel CGAL::Exact_predicates_inexact_constructions_kernel (recommended) 
\hasModel CGAL::Exact_predicates_exact_constructions_kernel 
\hasModel CGAL::Filtered_kernel 
\hasModel CGAL::Cartesian 
\hasModel CGAL::Simple_cartesian 
\hasModel CGAL::Homogeneous 
\hasModel CGAL::Simple_homogeneous 

*/

class FixedAlphaShapeTraits_3 {
public:

/// \name Types 
/// @{

/*! 
`CGAL::Comparison_result` or `Uncertain<CGAL::Comparison_result>` 
*/ 
typedef Hidden_type Comparison_result; 

/*! 
An object constructor able to compare the squared radius of the smallest circumscribing sphere of 
either four, three, two or one point(s) 
with a given value of alpha. 
It provides: 
- `Comparison_result operator()(Point_3,Point_3,Point_3,Point_3)` 
- `Comparison_result operator()(Point_3,Point_3,Point_3)` 
- `Comparison_result operator()(Point_3,Point_3)` 
- `Comparison_result operator()(Point_3)` 

*/ 
typedef Hidden_type Compare_squared_radius_3; 

/// @} 

/// \name Creation 
/// @{

/*! 
Default constructor. 
*/ 
FixedAlphaShapeTraits_3(); 

/// @} 

/// \name Access Functions 
/// @{

/*! 

*/ 
Compare_squared_radius_3 compare_squared_radius_3_object(); 

/// @}

}; /* end FixedAlphaShapeTraits_3 */

