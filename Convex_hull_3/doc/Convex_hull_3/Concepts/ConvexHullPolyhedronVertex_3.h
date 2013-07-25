/*!
\ingroup PkgConvexHull3Concepts
\cgalConcept

Requirements of the vertex type of the polyhedron built by the 
function `CGAL::convex_hull_3()`. 

\cgalHasModel CGAL::Polyhedron_3::Vertex

\sa `CGAL::Polyhedron_3` 
\sa `ConvexHullPolyhedronFacet_3` 
\sa `ConvexHullPolyhedronHalfedge_3` 

*/
class ConvexHullPolyhedronVertex_3 {
public:

/// \name Creation 
/// @{

/*!
default constructor. 
*/ 
Vertex(); 

/// @} 

/// \name Operations 
/// @{

/*!

*/ 
Point& point(); 

/*!
the point. 
*/ 
const Point& point() const; 

/// @}

}; /* end ConvexHullPolyhedronVertex_3 */
