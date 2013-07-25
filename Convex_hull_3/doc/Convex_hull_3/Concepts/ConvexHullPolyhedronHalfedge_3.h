/*!
\ingroup PkgConvexHull3Concepts
\cgalConcept

Requirements of the halfedge type required for the polyhedron built 
by the function `CGAL::convex_hull_3()`. 

\cgalHasModel CGAL::Polyhedron_3::Halfedge 

\sa `CGAL::Polyhedron_3` 
\sa `ConvexHullPolyhedronVertex_3` 
\sa `ConvexHullPolyhedronFacet_3` 

*/

class ConvexHullPolyhedronHalfedge_3 {
public:

/// \name Creation 
/// @{

/*!
default constructor. 
*/ 
Halfedge(); 

/// @} 

/// \name Operations 
/// @{

/*!
the opposite halfedge. 
*/ 
Halfedge_handle opposite(); 

/*!
the next halfedge around the facet. 
*/ 
Halfedge_handle next(); 

/*!
the previous halfedge around the facet. 
*/ 
Halfedge_handle prev(); 

/*!
is true if `h` is a border halfedge. 
*/ 
bool is_border() const; 

/*!
circulator of halfedges around the facet (counterclockwise). 
*/ 
Halfedge_around_facet_circulator facet_begin(); 

/*!
the incident vertex of `h`. 
*/ 
Vertex_handle vertex(); 

/*!
the incident facet of `h`. If `h` is a border halfedge 
the result is default construction of the handle. 
*/ 
Facet_handle facet(); 

/// @}

}; /* end ConvexHullPolyhedronHalfedge_3 */
