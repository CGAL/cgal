
namespace CGAL {
namespace Surface_mesh_simplification {

/*!
\ingroup PkgSurfaceMeshSimplificationRef

The class `Constrained_placement` is a model for the concept `GetPlacement`.
The placement of the vertex resulting from a contraction of an edge adjacent to a constrained edge
is the point of the common vertex. Otherwise the placement is the one computed by `Get_placement_`.

\tparam Get_placement_ must be a model of `GetPlacement`.
\tparam Edge_is_constrained_map_ must be a model of `ReadablePropertyMap` with `Edge_profile::edge_descriptor`
                                 as key type and `bool` as value type indicating if an edge is constrained.

\cgalModels `GetPlacement`

*/
template<class Get_placement_, class Edge_is_constrained_map_>
class Constrained_placement
  : public Get_placement_
{
public:

  /// \name Creation
  /// @{

  /*!
  Constructor
  */
  Constrained_placement(Edge_is_constrained_map_ map = Edge_is_constrained_map_(),
                        const Get_placement_& get_placement = Get_placement_());
  /// @}

  boost::optional<typename Edge_profile::Point> operator()(const Edge_profile& profile) const

};

} // namespace Surface_Mesh_Simplification
} // namespace CGAL
