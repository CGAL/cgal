namespace CGAL {
namespace Surface_mesh_simplification {

/*!
\ingroup PkgSurfaceMeshSimplificationRef

The class `GarlandHeckbert_plane_and_line_policies` regroups the cost and placement policies
based on the Garland-Heckbert "Plane and line" strategy of Liu and colleagues \cgalCite{liu2025linequadrics}.

This policy enhances the original Garland-Heckbert quadric error metrics
by adding to the cost the distance to the line passing through the input vertices and aligned with their normals.
Compared to the "classic" plane strategy, this strategy improves the speed and the quality of the result
(see Section \ref SurfaceMeshSimplificationGarlandHeckbertStrategy).

\note Both the cost and the placement policies must be used together as they internally use
and share information associating quadrics to vertices.
Note however, that they may still be wrapped with behavior-modifying classes
such as `Constrained_placement` or `Bounded_normal_change_placement`.

\tparam TriangleMesh is the type of surface mesh being simplified, and must be a model
                     of the `MutableFaceGraph` and `HalfedgeListGraph` concepts.
\tparam GeomTraits must be a model of `Kernel`. If you have passed a traits class in the optional
                   named parameters in the call to `edge_collapse()`, the types must be identical.

These policies depend on the third party \ref thirdpartyEigen library.

\sa `GarlandHeckbert_policies`
\sa `GarlandHeckbert_plane_policies`
\sa `GarlandHeckbert_probabilistic_plane_policies`
\sa `GarlandHeckbert_triangle_policies`
\sa `GarlandHeckbert_probabilistic_triangle_policies`
*/
template <typename TriangleMesh, typename GeomTraits>
class GarlandHeckbert_plane_and_line_policies
{
public:
  /// The type of the Garland-Heckbert cost functor, a model of the concept `GetCost`
  typedef unspecified_type Get_cost;

  /// The type of the Garland-Heckbert placement functor, a model of the concept `GetPlacement`
  typedef unspecified_type Get_placement;

  /// \name Creation
  /// @{

  /*!
  initializes the Garland-Heckbert plane and line policies.

  \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"

  \param tmesh the triangle mesh
  \param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below

  \cgalNamedParamsBegin
    \cgalParamNBegin{vertex_normal_map}
      \cgalParamDescription{a property map associating to each vertex of `tmesh` a normal direction for that vertex}
      \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<PolygonMesh>::%vertex_descriptor`
                     as key type and `Vector_3` as value type}
      \cgalParamDefault{an internal map filled using `CGAL::Polygon_mesh_processing::compute_vertex_normals`}
    \cgalParamNEnd
    \cgalParamNBegin{discontinuity_multiplier}
      \cgalParamDescription{a multiplier of the error value for boundary edges to preserve the boundaries}
      \cgalParamType{double}
      \cgalParamDefault{100}
    \cgalParamNEnd
    \cgalParamNBegin{line_policies_weight}
      \cgalParamDescription{a value that defines the weight of the line policies compared to the plane policies}
      \cgalParamType{double}
      \cgalParamDefault{0.01}
    \cgalParamNEnd
  \cgalNamedParamsEnd
  */
  template<typename NamedParameters = parameters::Default_named_parameters>
  GarlandHeckbert_plane_and_line_policies(TriangleMesh& tmesh, const NamedParameters &np = parameters::default_values());

  /// @}

  /// \name Accessors
  /// @{

  const Get_cost& get_cost() const;
  const Get_placement& get_placement() const;

  /// @}

};

/*!
creates a Garland-Heckbert plane and line policies object.
*/
template<typename TriangleMesh, typename NamedParameters = parameters::Default_named_parameters>
GarlandHeckbert_plane_and_line_policies make_GarlandHeckbert_plane_and_line_policies(TriangleMesh& tmesh,
                                                                                     const NamedParameters& np = parameters::default_values());

} // namespace Surface_mesh_simplification
} // namespace CGAL
