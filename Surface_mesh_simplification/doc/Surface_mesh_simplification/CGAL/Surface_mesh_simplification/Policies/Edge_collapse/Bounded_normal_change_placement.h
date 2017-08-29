
namespace CGAL {
namespace Surface_mesh_simplification {

/*!
\ingroup PkgSurfaceMeshSimplification

The class `Bounded_normal_change_placement` is a model for the `GetPlacement` concept
which serves as a filter for another placement. It rejects the placement if any
triangle in the profile changes the normal by more than 90 degree.

\tparam Placement must be a model of the concept `GetPlacement`. 

\cgalModels `GetPlacement`

*/
template< typename Placement >
class Bounded_normal_change_placement {
public:

/// \name Creation 
/// @{

/*!
Default constructor 
*/ 
Bounded_normal_change_placement<Placement>(); 

/*!
Constructor 

@param place is the placement that will be filtered.
*/ 
Bounded_normal_change_placement<Placement>(const Placement& place); 

/// @} 

/// \name Operations 
/// @{

/*!
Returns the placement computed by `place`, if no 
triangle in the profile changes the normal by more than 90 degree.
*/ 
template <typename Profile>
optional<typename Profile::Point> 
operator()( Profile const& profile ) const; 

/// @}

}; /* end Surface_mesh_simplification::Bounded_normal_change_placement */
} /* end namespace Surface_Mesh_Simplification */
} /* end namespace CGAL */
