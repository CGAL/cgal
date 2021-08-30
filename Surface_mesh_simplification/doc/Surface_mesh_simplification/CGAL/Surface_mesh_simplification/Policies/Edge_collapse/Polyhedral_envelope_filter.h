
namespace CGAL {
namespace Surface_mesh_simplification {

/*!
\ingroup PkgSurfaceMeshSimplificationRef

The class `Polyhedral_envelope_filter` is a model for the `PlacementFilter` concept.
It rejects the placement if the nested filter rejects it, or
if any triangle in the profile is not inside the polyhedral envelope, in this order.


\tparam Filter must be a model of the concept `PlacementFilter`.  It defaults to a class that does
not filter any placement.

\cgalModels `PlacementFilter`

\sa `Polyhedral_envelope`

*/
template <typename GeomTraits, typename Filter>
class Polyhedral_envelope_filter {
public:

  /// The number type
  typedef typename Geom_traits::FT    FT;

  /// \name Creation
  /// @{

  /*!
  %Default constructor
  */
  Polyhedral_envelope_filter();

  /*!
  Constructor

  \param dist is the parameter given to the polyhedral envelope
  \param filter is the filter that will be filtered.
  */
  Polyhedral_envelope_filter(const FT& dist, const Filter& filter);

  /// @}

  /// \name Operations
  /// @{

  /*!
  returns the placement, if it does not get filtered by the wrapped filter,
  and if all triangles in the profile are inside the polyhedral envelope.
  */
  boost::optional<typename Edge_profile::Point> operator()(const Edge_profile& profile,
                                                           boost::optional<typename Profile::Point> op) const;

  /// @}

}; /* end Surface_mesh_simplification::Polyhedral_envelope_filter */
} // namespace Surface_Mesh_Simplification
} // namespace CGAL
