
namespace CGAL {
namespace Surface_mesh_simplification {

/*!
\ingroup PkgSurfaceMeshSimplificationRef

The class `Constrained_placement` is a model for the `GetPlacement` concept
provided the template parameter `BasePlacement` is such a model.
The placement of the vertex resulting from a contraction of an edge adjacent to a constrained edge
is the point of the common vertex. Otherwise the placement is the one computed by `BasePlacement`.

\tparam BasePlacement a model of `GetPlacement`.
\tparam EdgeIsConstrainedMap a model of `ReadablePropertyMap` with `Edge_profile::edge_descriptor`
                             as key type and `bool` as value type indicating if an edge is constrained.

\cgalModels `GetPlacement`

*/
template<class BasePlacement, class EdgeIsConstrainedMap>
class Constrained_placement
  : public BasePlacement
{
public:

  /// \name Creation
  /// @{

  /*!
  Constructor
  */
  Constrained_placement(EdgeIsConstrainedMap map = EdgeIsConstrainedMap(),
                        BasePlacement base_placement = BasePlacement());
  /// @}

  boost::optional<typename Edge_profile::Point> operator()(const Edge_profile& profile) const

};

} // namespace Surface_Mesh_Simplification
} // namespace CGAL
