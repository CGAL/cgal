namespace CGAL {

/*!
\ingroup PkgMesh_3MeshClasses

The class `Mesh_triangulation_3` is a metafunctor which provides the triangulation type to be used
for the 3D triangulation embedding the mesh.

\tparam MD must be a model of `MeshDomain_3`.

\tparam Gt must be a model of `RegularTriangulationTraits_3` or `Default`
and defaults to `Kernel_traits<MD>::%Kernel`.

\tparam Concurrency_tag enables sequential versus parallel meshing and optimization algorithms.
                        Possible values are `Sequential_tag` (the default) and
                        `Parallel_tag`.

\tparam Vertex_base must be a model of `MeshVertexBase_3` or `Default`
and defaults to `Mesh_vertex_base_3<Gt, MD>`.

\tparam Cell_base must be a model of `MeshCellBase_3` or `Default`
and defaults to `Compact_mesh_cell_base_3<Gt, MD>`.

\warning To improve the robustness of the meshing process, the input traits `Gt`
         is wrapped with the traits class `Robust_weighted_circumcenter_filtered_traits_3`.
         The class `Robust_weighted_circumcenter_filtered_traits_3<Gt>` upgrades the functors
         models of `Kernel::ConstructWeightedCircumcenter_3`, `Kernel::ComputeSquaredRadius_3`,
         and `Kernel::ComputeSquaredRadiusSmallestOrthogonalSphere_3` that are
         provided by `Gt` to use exact computations when the geometric configuration
         is close to degenerate (e.g. almost coplanar points). <br>
         Users should therefore be aware that the traits class of the triangulation
         will have type `Robust_weighted_circumcenter_filtered_traits_3<Gt>`.

\sa `make_mesh_3()`
\sa `Mesh_complex_3_in_triangulation_3<Tr,CornerIndex,CurveSegmentIndex>`

*/
template< typename MD, typename Gt,
          typename Concurrency_tag,
          typename Vertex_base,
          typename Cell_base >
class Mesh_triangulation_3 {
public:

/// \name Types
/// @{

/*!
The triangulation type to be used
for the 3D triangulation embedding the mesh.
This type is a `Regular_triangulation_3` type
whose vertex and cell base classes are respectively
`Vertex_base` and `Cell_base`.
*/
typedef unspecified_type type;

/// @}

}; /* end Mesh_triangulation_3 */
} /* end namespace CGAL */
