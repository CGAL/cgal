
namespace CGAL {
namespace Surface_mesh_simplification {

/*!
\ingroup PkgSurfaceMeshSimplification

The class `Count_ratio_stop_predicate` is a model for the `StopPredicate` concept
which returns `true` when the relation between the initial and current number of edges drops below a certain ratio.  

\tparam ECM is the type of surface mesh being simplified, and must be a model of the `EdgeCollapsableSurfaceMesh` concept. 


\cgalModels `StopPredicate`

\sa `CGAL::Surface_mesh_simplification::Count_stop_predicate<ECM>` 

*/
template< typename ECM >
class Count_ratio_stop_predicate {
public:

/// \name Creation 
/// @{

/*!
Initializes the predicate establishing the `ratio`. 
*/ 
Count_ratio_stop_predicate<ECM>( double ratio ); 

/// @} 

/// \name Operations 
/// @{

/*!
Returns ` ( ((double)current_count / (double)initial_count) < ratio)`. 
All other parameters are ignored (but exist since this is a generic policy). 
*/ 
bool operator()( FT const& current_cost 
, Profile const& edge_profile 
, size_type initial_count 
, size_type current_count 
) const ; 

/// @}

}; /* end Surface_mesh_simplification::Count_ratio_stop_predicate */
} /* namespace Surface_mesh_simplification */
} /* end namespace CGAL */
