namespace CGAL {
namespace Surface_mesh_simplification {
/*!
\ingroup PkgSurfaceMeshSimplificationRef

The class `GarlandHeckbert_cost` provides a model for the `GetCost` concept.
It computes the collapse cost following the Garland-Heckbert strategy
(Section \ref SurfaceMeshSimplificationGarlandHeckbertStrategy).

It must be used in conjonction with the Garland-Heckbert placement policy,
`CGAL::Surface_mesh_simplification::GarlandHeckbert_placement<TriangleMesh>`

\tparam VertexCostMap must be a model of `ReadWritePropertyMap` with value type ``.

\cgalModels `GetCost`

\sa `CGAL::Surface_mesh_simplification::GarlandHeckbert_placement<TriangleMesh>`

*/
template<class VertexCostMap>
class GarlandHeckbert_cost
{
public:
  /// \name Creation
  /// @{

  /*!
  Initializes the policy with the given <I>garland heckbert state</I> object.
  Garland&Heckbert strategy requires a shared state object between cost and placementr policies.
  */
  GarlandHeckbert_cost(const VertexCostMap&);
  /// @}


  /// \name Operations
  /// @{

  /*!
  Returns the cost of collapsing the edge (represented by its profile) considering
  the new `placement` computed for it.
  */
  boost::optional<typename Edge_profile::FT>
  operator()(const Edge_profile& profile,
             const boost::optional<typename Edge_profile::Point>& placement) const;
  /// @}
};

} // namespace Surface_mesh_simplification
} // namespace CGAL
