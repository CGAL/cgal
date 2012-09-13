
namespace CGAL {

/*!
\ingroup PkgSurfaceMeshSimplification

The class `Surface_mesh_simplification::Midpoint_placement` provides a model for the `GetPlacement` concept. 
It computes the placement as the midpoint position along the edge. 

The class `Surface_mesh_simplification::Midpoint_placement` has one template arguments: the type of surface being simplified. 
It be a model of the `EdgeCollapsableMesh` concept. 

\models ::GetPlacement 

*/
template< typename ECM >
class Surface_mesh_simplification::Midpoint_placement {
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
result_type operator()( Profile const& edge_profile ) const; 

/// @}

}; /* end Surface_mesh_simplification::Midpoint_placement */
} /* end namespace CGAL */
