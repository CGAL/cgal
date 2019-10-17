
namespace CGAL {
namespace Surface_mesh_simplification {

/*!
\ingroup PkgSurfaceMeshSimplificationRef

The class `Bounded_normal_change_placement` is a model for the `GetPlacement` concept
which serves as a filter for another placement. It rejects the placement if any
triangle in the profile changes the normal by more than 90 degree.

\tparam Placement must be a model of the concept `GetPlacement`.

\cgalModels `GetPlacement`

*/
template <typename Placement>
class Bounded_normal_change_placement {
public:

/// \name Creation
/// @{

/*!
%Default constructor
*/
Bounded_normal_change_placement<Placement>();

/*!
Constructor

\param `base_placement` is the placement that will be filtered.
*/
Bounded_normal_change_placement<Placement>(const Placement& base_placement);

/// @}

/// \name Operations
/// @{

/*!
Returns the placement computed by `base_placement`, if no
triangle in the profile has its normal change by more than 90 degree.
*/
boost::optional<typename Edge_profile::Point> operator()(const Edge_profile& profile) const;

/// @}

}; /* end Surface_mesh_simplification::Bounded_normal_change_placement */
} // namespace Surface_Mesh_Simplification
} // namespace CGAL
