namespace CGAL {
namespace Surface_mesh_simplification {

/*!
\ingroup PkgSurfaceMeshSimplificationRef

This class is a type accessor for the quadric type.
*/
template <typename GeomTraits>
struct GarlandHeckbert_cost_matrix
{
  /// The type of a quadric in the Garland-Heckbert algorithm
  typedef typename Eigen::Matrix<typename GeomTraits::FT, 4, 4>           type;
};

/*!
\ingroup PkgSurfaceMeshSimplificationRef

The class `GarlandHeckbert_cost` provides a model for the `GetCost` concept.
It computes the collapse cost following the Garland-Heckbert strategy
(Section \ref SurfaceMeshSimplificationGarlandHeckbertStrategy).

It must be used in conjonction with the Garland-Heckbert placement policy,
`CGAL::Surface_mesh_simplification::GarlandHeckbert_placement<TriangleMesh>`

\tparam VertexCostMap must be a model of `ReadWritePropertyMap` with value type `GarlandHeckbert_cost_matrix`.

\cgalModels `GetCost`

\sa `CGAL::Surface_mesh_simplification::GarlandHeckbert_placement<VertexCostMap>`

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
