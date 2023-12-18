
namespace CGAL {
namespace Surface_mesh_simplification {

/*!
\ingroup PkgSurfaceMeshSimplificationRef

\deprecated This class is deprecated since \cgal 5.3 and the use of
`Bounded_normal_change_filter` should be preferred.


The class `Bounded_normal_change_placement` is a model for the `GetPlacement` concept
which serves as a filter for another placement. It rejects the placement if any
triangle in the profile changes the normal by more than 90 degree.

\tparam Get_placement_ must be a model of the concept `GetPlacement`.

\cgalModels{GetPlacement}

*/
template <typename Get_placement_>
class Bounded_normal_change_placement {
public:

  /// \name Creation
  /// @{

  /*!
  %Default constructor
  */
  Bounded_normal_change_placement();

  /*!
  Constructor

  \param get_placement is the placement that will be filtered.
  */
  Bounded_normal_change_placement(const Get_placement_& get_placement);

  /// @}

  /// \name Operations
  /// @{

  /*!
  Returns the placement computed by `get_placement`, if no
  triangle in the profile has its normal changed by more than 90 degree.
  */
  std::optional<typename Edge_profile::Point> operator()(const Edge_profile& profile) const;

  /// @}

}; /* end Surface_mesh_simplification::Bounded_normal_change_placement */
} // namespace Surface_Mesh_Simplification
} // namespace CGAL
