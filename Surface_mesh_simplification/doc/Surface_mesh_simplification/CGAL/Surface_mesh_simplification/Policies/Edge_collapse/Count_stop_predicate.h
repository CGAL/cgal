
namespace CGAL {
namespace Surface_mesh_simplification {

/*!
\ingroup PkgSurfaceMeshSimplification

The class `Count_stop_predicate` is a model for the `StopPredicate` concept,
which returns `true` when the number of current edges drops below a certain threshold. 
 
\tparam ECM is the type of surface mesh being simplified, and must be a model of the `EdgeCollapsableSurfaceMesh` concept. 

\cgalModels `StopPredicate`

\sa `CGAL::Surface_mesh_simplification::Count_ratio_stop_predicate<ECM>` 

*/
template< typename ECM >
class Count_stop_predicate {
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
Returns `(current_count < threshold)`. All other parameters are ignored (but exist since this is a generic policy). 
*/ 
bool operator()( FT const& current_cost 
, Profile const& edge_profile 
, size_type initial_count 
, size_type current_count 
) const ; 

/// @}

}; /* end Surface_mesh_simplification::Count_stop_predicate */
} /* namespace Surface_mesh_simplification */
} /* end namespace CGAL */
