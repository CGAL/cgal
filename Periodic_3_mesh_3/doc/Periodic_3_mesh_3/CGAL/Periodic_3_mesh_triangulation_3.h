namespace CGAL {

/*!
\ingroup PkgPeriodic_3_mesh_3MeshClasses

The class `Periodic_3_mesh_triangulation_3` is a class template which provides
the triangulation type to be used for the 3D periodic triangulation embedding the mesh.

\tparam MD must be a model of `Periodic_3MeshDomain_3` or `Periodic_3MeshDomainWithFeatures_3`.

\tparam Gt must be a model of `MeshTriangulationTraits_3`.
`Default` may be used, with default value `Kernel_traits<MD>::%Kernel`.

\tparam Vertex_base must be a model of `MeshVertexBase_3` and `Periodic_3TriangulationDSVertexBase_3`.
`Default` may be used, with default value `Mesh_vertex_base_3<Gt, MD, Triangulation_vertex_base_3<Gt, Periodic_3_triangulation_ds_vertex_base_3> >`.

\tparam Cell_base must be a model of `MeshCellBase_3` and `Periodic_3TriangulationDSCellBase_3`.
`Default` may be used, with default value `Mesh_cell_base_3<Gt, MD, Triangulation_cell_base_with_circumcenter_3<Gt, Triangulation_cell_base_3<Gt, Periodic_3_triangulation_ds_cell_base_3> > >`.

\warning The input traits `Gt` are wrapped multiple times to handle periodicity
         and to improve the robustness of the meshing process: for example,
         wrapping `Gt` with `Robust_weighted_circumcenter_filtered_traits_3<Gt>`
         allows the functors models of `Kernel::ConstructWeightedCircumcenter_3`,
         `Kernel::ComputeSquaredRadius_3`, and `Kernel::ComputeSquaredRadiusSmallestOrthogonalSphere_3`
          that are provided by `Gt` to use exact computations when the geometric configuration
         is close to degenerate (e.g. almost coplanar points). <br>
         Users should therefore be aware that the type of the traits class
         of the triangulation will not be `Gt`.

\sa `make_periodic_3_mesh_3()`
\sa `refine_periodic_3_mesh_3()`

\sa `Mesh_triangulation_3`

*/
template< typename MD, typename Gt,
          typename Vertex_base,
          typename Cell_base>
class Periodic_3_mesh_triangulation_3
{
public:

  /// \name Types
  /// @{

  /*!
  The triangulation type to be used for the 3D triangulation embedding the mesh.
  This type is a wrapper around the type `CGAL::Periodic_3_regular_triangulation_3`,
  whose vertex and cell base classes are respectively `Vertex_base` and `Cell_base`.
  */
  typedef unspecified_type type;

  /// @}

}; /* end Periodic_3_mesh_triangulation_3 */
} /* end namespace CGAL */
