
namespace CGAL {

/*!
\ingroup PkgSurfaceMeshSimplification

The class `Surface_mesh_simplification::Count_ratio_stop_predicate` provides a model for the `StopPredicate` concept. 
It has one template argument: the type of surface being simplified, 
which must be a model of the `EdgeCollapsableMesh` concept. 
It returns `true` when the relation between the initial and current number 
of edges drops below a certain ratio. 

\models ::StopPredicate 

\sa `CGAL::Surface_mesh_simplification::Count_stop_predicate<ECM>` 

*/
template< typename ECM >
class Surface_mesh_simplification::Count_ratio_stop_predicate {
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
Returns \f$ ( ((double)current\_count / (double)initial\_count) < ratio)\f$. 
All other parameters are ignored (but exist since this is a generic policy). 
*/ 
bool operator()( FT const& current_cost 
, Profile const& edge_profile 
, size_type initial_count 
, size_type current_count 
) const ; 

/// @}

}; /* end Surface_mesh_simplification::Count_ratio_stop_predicate */
} /* end namespace CGAL */
