namespace CGAL {
namespace Surface_mesh_simplification {
/*!
\ingroup PkgSurfaceMeshSimplificationRef

The class `GarlandHeckbert_cost` provides a model for the `GetCost` concept.
It computes the collapse cost following the Garland-Heckbert strategy
(Section \ref SurfaceMeshSimplificationGarlandHeckbertStrategy)

\tparam TriangleMesh is the type of surface mesh being simplified, and must be a model of the `MutableFaceGraph` and `HalfedgeListGraph` concepts.

\cgalModels `GetCost`

\sa `CGAL::Surface_mesh_simplification::GarlandHeckbert_placement<TriangleMesh>`

*/
template<class TriangleMesh>
class GarlandHeckbert_cost
{
public:
  /// \name Creation
  /// @{

  /*!
  Initializes the policy with the given <I>garland heckbert state</I> object.
  Garland&Heckbert strategy requires a shared state object between cost, placement, and visitor policies.
  */
  GarlandHeckbert_cost();

  /// @}


  /// \name Operations
  /// @{

  /*!
  Returns the cost of collapsing the edge (represented by its profile) considering
  the new `placement` computed for it.
  */
  template <typename Profile>
  boost::optional<typename Profile::FT>
  operator()(const Profile& profile,
             const boost::optional<typename Profile::Point>& placement) const;
  /// @}
};/* end Surface_mesh_simplification::GarlandHeckbert_cost */

} // namespace Surface_mesh_simplification
} // namespace CGAL
