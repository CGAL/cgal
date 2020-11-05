namespace CGAL {
namespace Surface_mesh_simplification {

/*!
\ingroup PkgSurfaceMeshSimplificationRef

The class `LindstromTurk_placement` provides a model for the `GetPlacement` concept.
It computes the placement, that is, the new position for the remaining vertex after
a halfedge collapse, following the Lindstrom-Turk strategy
(Section \ref SurfaceMeshSimplificationLindstromTurkStrategy).

\tparam TriangleMesh is the type of surface mesh being simplified, and must be a model of the `MutableFaceGraph` and `HalfedgeListGraph` concepts.

\cgalModels `GetPlacement`

\sa `CGAL::Surface_mesh_simplification::LindstromTurk_cost<TriangleMesh>`

*/
template <typename TriangleMesh>
class LindstromTurk_placement
{
public:
  /// \name Creation
  /// @{

  /*!
  Initializes the policy with the given <I>weighting unit factor</I>.
  See \ref SurfaceMeshSimplificationLindstromTurkStrategy for details on the meaning of this factor.
  */
  LindstromTurk_placement(const Edge_profile::FT& factor = FT(0.5));

  /// @}

  /// \name Operations
  /// @{

  /*!
  Returns the new position for the remaining vertex after collapsing the edge
  (represented by its profile).
  */
  boost::optional<typename Edge_profile::Point>
  operator()(const Edge_profile& profile) const;

  /// @}

};

} // namespace Surface_mesh_simplification
} // namespace CGAL
