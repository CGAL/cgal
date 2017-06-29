
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

\cgalRefines `CopyConstructible` 

\cgalHasModel `CGAL::Surface_mesh_simplification::Midpoint_placement<ECM>`
\cgalHasModel `CGAL::Surface_mesh_simplification::LindstromTurk_placement<ECM>`
\cgalHasModel `CGAL::Surface_mesh_simplification::Bounded_normal_change_placement<Placement>`

*/

class GetPlacement {
public:


/// \name Operations 
/// @{

/*!
Computes and returns the placement, that is, the position of the vertex 
which replaces the collapsing edge (represented by its profile). 

\tparam Profile must be a model of `EdgeProfile`.
*/ 
template <class Profile>
boost::optional<Point>operator()( Profile const& edge_profile ) const; 

/// @}

}; /* end GetPlacement */

