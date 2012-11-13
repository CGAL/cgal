
/*!
\ingroup PkgSurfaceMeshSimplificationConcepts
\cgalConcept

The concept `GetPlacement` describes the requirements for the <I>policy 
function object</I> which gets the <I>collapse placement</I> of an edge, 
that is, the new position of the vertex that remains after a 
halfedge-collapse operation. 

The placement returned is a `boost::optional` value (i.e., it can 
be absent). An absent result indicates that the remaining vertex 
must be kept in place, not moved to a new position. 

\cgalRefines `DefaultConstructible` 
\cgalRefines `CopyConstructible` 

\cgalHasModel `CGAL::Surface_mesh_simplification::Midpoint_placement<ECM>`
\cgalHasModel `CGAL::Surface_mesh_simplification::LindstromTurk_placement<ECM>`

*/

class GetPlacement {
public:

/// \name Types 
/// @{

/*! 
The type of the edge profile cache. Must be a model of the `EdgeProfile` concept. 
*/ 
typedef Hidden_type Profile; 

/*! 
The point type for the surface vertex. Must be a model of `Point_3`. 
*/ 
typename CGAL::halfedge_graph_traits<ECM>::Point Point; 

/*! 
The type of the result (an optional point). 
*/ 
boost::optional<Point> result_type; 

/// @} 

/// \name Operations 
/// @{

/*! 
Computes and returns the placement, that is, the position of the vertex 
which replaces the collapsing edge (represented by its profile). 
*/ 
result_type operator()( Profile const& edge_profile ) const; 

/// @}

}; /* end GetPlacement */

