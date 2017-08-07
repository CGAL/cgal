namespace CGAL {

/*!
\ingroup PkgPeriodic_3_mesh_3MeshClasses

The class `Periodic_3_mesh_triangulation_3` is a metafunctor which provides
the triangulation type to be used for the 3D periodic triangulation embedding the mesh.

\tparam MD must be a model of `Periodic_3MeshDomain_3`.

\tparam Gt must be a model of `RegularTriangulationTraits_3`.
`Default` may be used, with default value `Kernel_traits<MD>::%Kernel`.

\tparam Vertex_base must be a model of `MeshVertexBase_3` and `Periodic_3TriangulationDSVertexBase_3`.
`Default` may be used, with default value `Mesh_vertex_base_3<Gt, MD, Triangulation_vertex_base_3<Gt, Periodic_3_triangulation_ds_vertex_base_3> >`.

\tparam Cell_base must be a model of `MeshCellBase_3` and `Periodic_3TriangulationDSCellBase_3`.
`Default` may be used, with default value `Mesh_cell_base_3<Gt, MD, Triangulation_cell_base_with_circumcenter_3<Gt, Triangulation_cell_base_3<Gt, Periodic_3_triangulation_ds_cell_base_3> > >`.

\sa `make_periodic_3_mesh_3()`
\sa `refine_periodic_3_mesh_3()`

\sa `Mesh_triangulation_3`

*/
template< typename MD, typename Gt,
          typename Vertex_base,
          typename Cell_base>
class Periodic_3_mesh_triangulation_3 {
public:

/// \name Types
/// @{

/*!
The triangulation type to be used for the 3D triangulation embedding the mesh.
This type is a `Periodic_3Regular_triangulation_3` type whose vertex and cell
base classes are respectively `Vertex_base` and `Cell_base`.
*/
typedef unspecified_type type;

/// @}

}; /* end Periodic_3_mesh_triangulation_3 */
} /* end namespace CGAL */
