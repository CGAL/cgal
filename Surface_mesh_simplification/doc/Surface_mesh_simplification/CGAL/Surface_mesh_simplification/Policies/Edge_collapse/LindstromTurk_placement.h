
namespace CGAL {
namespace Surface_mesh_simplification {

/*!
\ingroup PkgSurfaceMeshSimplification

The class `LindstromTurk_placement` provides a model for the `GetPlacement` concept. 
It computes the placement, that is, the new position for the remaining vertex after 
a halfedge-collapse, following the Lindstrom-Turk strategy 
(Section \ref SurfaceMeshSimplificationLindstromTurkStrategy). 

\tparam ECM is the type of surface mesh being simplified, and must be a model of the `EdgeCollapsableSurfaceMesh` concept. 

\cgalModels `GetPlacement`

\sa `CGAL::Surface_mesh_simplification::LindstromTurk_cost<ECM>` 

*/
template< typename ECM >
class LindstromTurk_placement {
public:

/// \name Creation 
/// @{

/*!
Initializes the policy with the given <I>weighting unit factor</I>. 
See \ref SurfaceMeshSimplificationLindstromTurkStrategy for details on the meaning of this factor. 
*/ 
LindstromTurk_placement<ECM>( FT const& factor = FT(0.5) ); 

/// @} 

/// \name Operations 
/// @{

/*!
Returns the new position for the remaining vertex after collapsing the edge 
(represented by its profile). 
*/   
template <typename Profile> 
optional<typename Profile::Point>
operator()( Profile const& profile ) const; 

/// @}

}; /* end Surface_mesh_simplification::LindstromTurk_placement */
} /* namespace Surface_mesh_simplification */
} /* end namespace CGAL */
