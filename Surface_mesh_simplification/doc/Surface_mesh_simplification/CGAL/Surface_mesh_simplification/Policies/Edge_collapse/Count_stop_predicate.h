namespace CGAL {
namespace Surface_mesh_simplification {

/*!
\ingroup PkgSurfaceMeshSimplificationRef

The class `Count_stop_predicate` is a model for the `StopPredicate` concept,
which returns `true` when the number of current edges drops below a certain threshold.

\tparam TriangleMesh is the type of surface mesh being simplified, and must be a model of the `MutableFaceGraph` and `HalfedgeListGraph` concepts.

\cgalModels `StopPredicate`

\sa `CGAL::Surface_mesh_simplification::Count_ratio_stop_predicate<TriangleMesh>`

*/
template< typename TriangleMesh >
class Count_stop_predicate {
public:

/// \name Creation
/// @{

/*!
Initializes the predicate establishing the `threshold` value.
*/
Count_stop_predicate<TriangleMesh>(size_type threshold);

/// @}

/// \name Operations
/// @{

/*!
Returns `(current_count < threshold)`. All other parameters are ignored (but exist since this is a generic policy).
*/
bool operator()(const FT& current_cost,
                const Profile& edge_profile,
                size_type initial_count,
                size_type current_count) const;

/// @}

}; /* end Surface_mesh_simplification::Count_stop_predicate */
} /* namespace Surface_mesh_simplification */
} /* end namespace CGAL */
