
namespace CGAL {

/*!
\ingroup PkgSurfaceMeshSimplification

The class `Surface_mesh_simplification::Count_stop_predicate` provides a model for the `StopPredicate` concept. 
It has one template argument: the type of surface being simplified, 
which must be a model of the `EdgeCollapsableMesh` concept. 
It returns `true` when the number of current edges drops below a certain threshold. 

\models ::StopPredicate 

\sa `CGAL::Surface_mesh_simplification::Count_ratio_stop_predicate<ECM>` 

*/
template< typename ECM >
class Surface_mesh_simplification::Count_stop_predicate {
public:

/// \name Creation 
/// @{

/*! 
Initializes the predicate establishing the `threshold` value. 
*/ 
Count_stop_predicate<ECM>( size_type threshold ); 

/// @} 

/// \name Operations 
/// @{

/*! 
Returns \f$ (current\_count < threshold)\f$. All other parameters are ignored (but exist since this is a generic policy). 
*/ 
bool operator()( FT const& current_cost 
, Profile const& edge_profile 
, size_type initial_count 
, size_type current_count 
) const ; 

/// @}

}; /* end Surface_mesh_simplification::Count_stop_predicate */
} /* end namespace CGAL */
