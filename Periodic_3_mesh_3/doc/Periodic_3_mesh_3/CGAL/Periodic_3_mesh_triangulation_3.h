namespace CGAL {

/*!
\ingroup PkgPeriodic3Mesh3MeshClasses

The class `Periodic_3_mesh_triangulation_3` is a class template which provides
the triangulation type to be used for the 3D periodic triangulation embedding the mesh.

\tparam MD must be a model of `Periodic_3MeshDomain_3` or `Periodic_3MeshDomainWithFeatures_3`.

\tparam GT must be a model of `MeshTriangulationTraits_3`.
`Default` may be used, with default value `Kernel_traits<MD>::%Kernel`.

\tparam Vertex_base must be a model of `MeshVertexBase_3` and `Periodic_3TriangulationDSVertexBase_3`.
`Default` may be used, with default value `Mesh_vertex_base_3<GT, MD, Regular_triangulation_vertex_base_3<Gt, Periodic_3_triangulation_ds_vertex_base_3> >`.

\tparam Cell_base must be a model of `MeshCellBase_3` and `Periodic_3TriangulationDSCellBase_3`.
`Default` may be used, with default value `Mesh_cell_base_3<GT, MD, Regular_triangulation_cell_base_with_weighted_circumcenter_3<Gt, Regular_triangulation_cell_base_3<Gt, Periodic_3_triangulation_ds_cell_base_3> > >`.

\warning The input traits `GT` are wrapped multiple times to handle periodicity
         and to improve the robustness of the meshing process: for example,
         wrapping `GT` with `Robust_weighted_circumcenter_filtered_traits_3<GT>`
         allows the functors models of `Kernel::ConstructWeightedCircumcenter_3`,
         `Kernel::ComputeSquaredRadius_3`, and `Kernel::ComputeSquaredRadiusSmallestOrthogonalSphere_3`
          that are provided by `GT` to use exact computations when the geometric configuration
         is close to degenerate (e.g. almost coplanar points). <br>
         Users should therefore be aware that the type of the traits class
         of the triangulation will not be `GT`.

\sa `make_periodic_3_mesh_3()`
\sa `refine_periodic_3_mesh_3()`

\sa `Mesh_triangulation_3`

*/
template< typename MD, typename GT, typename Vertex_base, typename Cell_base>
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
