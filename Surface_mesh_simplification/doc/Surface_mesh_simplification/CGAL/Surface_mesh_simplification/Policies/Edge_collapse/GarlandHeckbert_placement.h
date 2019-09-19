namespace CGAL {
namespace Surface_mesh_simplification {

/*!
\ingroup PkgSurfaceMeshSimplificationRef

The class `GarlandHeckbert_placement` provides a model for the `GetPlacement` concept.
It computes the placement, that is, the new position for the remaining vertex after
a halfedge-collapse, following the Garland-Heckbert strategy
(Section \ref SurfaceMeshSimplificationGarlandHeckbertStrategy).

\tparam TriangleMesh is the type of surface mesh being simplified, and must be a model of the `MutableFaceGraph` and `HalfedgeListGraph` concepts.

\cgalModels `GetPlacement`

\sa `CGAL::Surface_mesh_simplification::GarlandHeckbert_cost<TriangleMesh>`

*/
template<class TriangleMesh>
class GarlandHeckbert_placement
{
public:
  /// \name Creation
  /// @{

  /*!
  Initializes the policy with the given <I>garland heckbert state</I> object.
  Garland&Heckbert strategy requires a shared state object between cost, placement, and visitor policies.
  */
  GarlandHeckbert_placement();
  /// @}

  /// \name Operations
  /// @{

  /*!
  Returns the new position for the remaining vertex after collapsing the edge
  (represented by its profile).
  */
  template <typename Profile>
  boost::optional<typename Profile::Point> operator()(const Profile& profile) const;
  /// @}
};
} // namespace Surface_mesh_simplification
} // namespace CGAL
