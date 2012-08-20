/*!
\ingroup PkgConvexHull3Concepts
\cgalconcept

The requirements of the vertex type of the polyhedron to be built by the 
function `::convex_hull_3`. 

\hasModel CGAL::Polyhedron_3<Traits>::Vertex

\sa `CGAL::Polyhedron_3<Traits>` 
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
