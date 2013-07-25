
/*!
\ingroup PkgAlphaShapes3Concepts
\cgalConcept

The concept `AlphaShapeTraits_3` describes the requirements 
for the geometric traits class 
of the underlying Delaunay triangulation of a basic alpha shape. 

\cgalRefines `DelaunayTriangulationTraits_3` 

In addition to the requirements described in the concept
::DelaunayTriangulationTraits_3, the geometric traits class of a
Delaunay triangulation plugged in a basic alpha shapes provides the
following.

\cgalHasModel All %CGAL kernels. 
\cgalHasModel `CGAL::Exact_predicates_inexact_constructions_kernel` (recommended) 
\cgalHasModel `CGAL::Exact_predicates_exact_constructions_kernel` 
\cgalHasModel `CGAL::Filtered_kernel` 
\cgalHasModel `CGAL::Cartesian` 
\cgalHasModel `CGAL::Simple_cartesian` 
\cgalHasModel `CGAL::Homogeneous` 
\cgalHasModel `CGAL::Simple_homogeneous`

*/

class AlphaShapeTraits_3 {
public:

/// \name Types 
/// @{

/*!
A number type compatible with the type used for 
the points coordinate. 
*/ 
typedef unspecified_type FT; 

/*!
An object constructor able to compute: 

- the squared radius of the 
smallest circumscribing sphere of 4 points `p0, p1, p2, p3`, 
- the squared radius of the 
smallest circumscribing sphere of 3 points `p0, p1, p2`, 
- the squared radius of the smallest circumscribing sphere 
of 2 points `p0, p1`, 
- the squared radius of the smallest circumscribing sphere 
to a single point `p0` (which is `FT(0)`). 
*/ 
typedef unspecified_type Compute_squared_radius_3; 

/// @} 

/// \name Creation 
/// @{

/*!
Default constructor. 
*/ 
AlphaShapeTraits_3(); 

/// @} 

/// \name Access Functions 
/// @{

/*!

*/ 
Compute_squared_radius_3 compute_squared_radius_3_object(); 

/// @}

}; /* end AlphaShapeTraits_3 */

