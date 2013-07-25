/*!
\ingroup PkgConvexHull3Concepts
\cgalConcept

Requirements of the facet type of a polyhedron built by the 
function `CGAL::convex_hull_3()`. 

\cgalHasModel `CGAL::Polyhedron_3::Facet` 

\sa `CGAL::Polyhedron_3` 
\sa `ConvexHullPolyhedronVertex_3` 
\sa `ConvexHullPolyhedronHalfedge_3` 

*/

class ConvexHullPolyhedronFacet_3 {
public:

/// \name Types 
/// @{

/*!
plane equation type stored in facets. 
*/ 
typedef unspecified_type Plane; 

/*!
handle to halfedge. 
*/ 
typedef unspecified_type Halfedge_handle; 

/*!
circulator of 
halfedges around a facet. 
*/ 
typedef unspecified_type Halfedge_around_facet_circulator; 

/// @} 

/// \name Creation 
/// @{

/*!
default constructor. 
*/ 
Facet(); 

/// @} 

/// \name Operations 
/// @{

/*!

*/ 
Plane& plane(); 

/*!
plane equation. 
*/ 
const Plane& plane() const; 

/*!
an incident halfedge that points to `f`. 
*/ 
Halfedge_handle halfedge(); 

/*!
circulator of halfedges around the facet (counterclockwise). 
*/ 
Halfedge_around_facet_circulator facet_begin(); 

/// @}

}; /* end ConvexHullPolyhedronFacet_3 */
