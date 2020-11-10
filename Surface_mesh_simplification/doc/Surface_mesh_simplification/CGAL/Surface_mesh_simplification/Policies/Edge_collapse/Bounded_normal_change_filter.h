
namespace CGAL {
namespace Surface_mesh_simplification {

/*!
\ingroup PkgSurfaceMeshSimplificationRef

The class `Bounded_normal_change_filter` is a model for the `Filter` concept.
It rejects the placement if the nested filter rejects it, or
if any triangle in the profile changes the normal by more than 90 degree, in this order.


\tparam Filter_ must be a model of the concept `Filter`.  It defaults to a class that does
not filter any placement.

\cgalModels `Filter`

*/
template <typename Filter_>
class Bounded_normal_change_filter {
public:

  /// \name Creation
  /// @{

  /*!
  %Default constructor
  */
  Bounded_normal_change_filter();

  /*!
  Constructor

  \param filter is the filter that will be filtered.
  */
  Bounded_normal_change_filter(const Filter_& filter);

  /// @}

  /// \name Operations
  /// @{

  /*!
  Returns the placement, if it does not get filtered by the wrapped filter
  and if no triangle in the profile has its normal changed by more than 90 degree.
  */
  boost::optional<typename Edge_profile::Point> operator()(const Edge_profile& profile,
                                                           boost::optional<typename Profile::Point> op) const;

  /// @}

}; /* end Surface_mesh_simplification::Bounded_normal_change_filter */
} // namespace Surface_Mesh_Simplification
} // namespace CGAL
