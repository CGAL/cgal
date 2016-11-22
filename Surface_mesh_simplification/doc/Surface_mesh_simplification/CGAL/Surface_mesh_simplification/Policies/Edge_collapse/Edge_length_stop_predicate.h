
namespace CGAL {
namespace Surface_mesh_simplification {

/*!
\ingroup PkgSurfaceMeshSimplification

The class `Edge_length_stop_predicate` is a model for the `StopPredicate` concept,
which returns `true` when the top edge in the priority queue is larger than a certain threshold.
This predicate is meant to be used with `Edge_length_cost`.

\tparam FT is the number type of the point coordinates.

\cgalModels `StopPredicate`

*/
template< typename FT >
class Edge_length_stop_predicate {
public:

/// \name Creation
/// @{

/*!
Initializes the predicate establishing the `threshold` value.
*/
Edge_length_stop_predicate<ECM>( FT threshold );

/// @}

/// \name Operations
/// @{

/*!
Returns `(CGAL::squared_distance(edge_profile.p0(),edge_profile.p1()) > threshold*threshold)`.
All other parameters are ignored (but exist since this is a generic policy).
*/
bool operator()( FT const&
, Profile const& edge_profile
, size_type
, size_type
) const ;

/// @}

}; /* end Surface_mesh_simplification::Edge_length_stop_predicate*/
} /* namespace Surface_mesh_simplification */
} /* end namespace CGAL */
