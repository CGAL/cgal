namespace CGAL {
namespace Surface_mesh_simplification {

/*!
\ingroup PkgSurfaceMeshSimplificationRef

The class `GarlandHeckbert_plane_policies` regroups the cost and placement policies
based on the Garland-Heckbert "Classic Plane" strategy
(Section \ref SurfaceMeshSimplificationGarlandHeckbertStrategy).
This class implements the original Garland-Heckbert quadric error metric strategy,
as described in their seminal paper \cgalCite{gh-ssqem-97}.

Both the cost and the placement policies must be used together as they internally use
and share information associating quadrics to vertices.
Note however, that they may still be wrapped with behavior modifying classes
such as `Constrained_placement` or `Bounded_normal_change_placement`.

\tparam TriangleMesh is the type of surface mesh being simplified, and must be a model
                     of the `MutableFaceGraph` and `HalfedgeListGraph` concepts.
\tparam GeomTraits must be a model of `Kernel`. If you have passed a traits class in the optional
                   named parameters in the call to `edge_collapse()`, the types must be identical.

These policies depend on the third party \ref thirdpartyEigen library.

\sa `GarlandHeckbert_probabilistic_plane_policies`
\sa `GarlandHeckbert_triangle_policies`
\sa `GarlandHeckbert_probabilistic_triangle_policies`
*/
template <typename TriangleMesh, typename GeomTraits>
class GarlandHeckbert_plane_policies
{
public:
  /// The type of the Garland Heckbert cost functor, a model of the concept `GetCost`
  typedef unspecified_type Get_cost;

  /// The type of the Garland Heckbert placement functor, a model of the concept `GetPlacement`
  typedef unspecified_type Get_placement;

  /// \name Creation
  /// @{

  /*!
  initializes the Garland-Heckbert Plane policies.
  */
  GarlandHeckbert_plane_policies(TriangleMesh& tmesh);

  /// @}

  /// \name Accessors
  /// @{

  const Get_cost& get_cost() const;
  const Get_placement& get_placement() const;

  /// @}

};

} // namespace Surface_mesh_simplification
} // namespace CGAL
