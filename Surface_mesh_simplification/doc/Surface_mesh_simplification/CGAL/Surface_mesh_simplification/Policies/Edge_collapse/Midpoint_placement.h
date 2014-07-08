
namespace CGAL {
namespace Surface_mesh_simplification {

/*!
\ingroup PkgSurfaceMeshSimplification

The class `Midpoint_placement` is a model for the `GetPlacement` concept
which  computes the placement as the midpoint position along the edge. 

\tparam ECM is the type of surface mesh being simplified, and must be a model of the `EdgeCollapsableSurfaceMesh` concept. 

\cgalModels `GetPlacement`

*/
template< typename ECM >
class Midpoint_placement {
public:

/// \name Creation 
/// @{

/*!
Default constructor 
*/ 
Midpoint_placement<ECM>(); 

/// @} 

/// \name Operations 
/// @{

/*!
Returns the <I>placement</I> (vertex position) as the midpoint between 
the points of the source and target vertices 
(`profile.p0()` and `profile.p1()`) 
*/ 
template <typename Profile>
optional<typename Profile::Point> 
operator()( Profile const& profile ) const; 

/// @}

}; /* end Surface_mesh_simplification::Midpoint_placement */
} /* end namespace Surface_Mesh_Simplification */
} /* end namespace CGAL */
