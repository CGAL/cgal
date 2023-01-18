namespace CGAL {
namespace Surface_mesh_simplification {

/*!
\ingroup PkgSurfaceMeshSimplificationRef

The class `GarlandHeckbert_probabilistic_plane_policies` regroups the cost and placement policies
based on the "Probabilistic Plane" strategy of Trettner and Kobbelt \cgalCite{tk-frqmpq-20}.
This policy enhances the original Garland-Heckbert quadric error metrics,
by adding notions of variance (i.e., noise) to the positions of the vertices and the normals of faces.
Compared to the "classic" plane strategy, this strategy sacrifices tightness to the input
for speed and quality of the result (Section \ref SurfaceMeshSimplificationGarlandHeckbertStrategy).

Both the cost and the placement policies must be used together as they internally use
and share information associating quadrics to vertices.
Note however, that they may still be wrapped with behavior modifying classes
such as `Constrained_placement` or `Bounded_normal_change_placement`.

\tparam TriangleMesh is the type of surface mesh being simplified, and must be a model
                     of the `MutableFaceGraph` and `HalfedgeListGraph` concepts.
\tparam GeomTraits must be a model of `Kernel`. If you have passed a traits class in the optional
                   named parameters in the call to `edge_collapse()`, the types must be identical.

These policies depend on the third party \ref thirdpartyEigen library.

\sa `GarlandHeckbert_plane_policies`
\sa `GarlandHeckbert_triangle_policies`
\sa `GarlandHeckbert_probabilistic_triangle_policies`
*/
template <typename TriangleMesh, typename GeomTraits>
class GarlandHeckbert_probabilistic_plane_policies
{
public:
  /// The type of the Garland Heckbert cost functor, a model of the concept `GetCost`
  typedef unspecified_type Get_cost;

  /// The type of the Garland Heckbert placement functor, a model of the concept `GetPlacement`
  typedef unspecified_type Get_placement;

  /// \name Creation
  /// @{

  /*!
  initializes the Garland-Heckbert Probabilistic Plane policies.
  */
  GarlandHeckbert_probabilistic_plane_policies(TriangleMesh& tmesh);

  /// @}

  /// \name Accessors
  /// @{

  const Get_cost& get_cost() const;
  const Get_placement& get_placement() const;

  /// @}
};

} // namespace Surface_mesh_simplification
} // namespace CGAL
