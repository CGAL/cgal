namespace CGAL {

/*!
\ingroup PkgMesh_3MeshClasses

The class `Compact_mesh_cell_base_3<Gt, MD>` is a model of the concept `MeshCellBase_3`.
It is designed to serve as cell base class for the 3D triangulation 
used in the 3D mesh generation process. It is more compact in memory than
`Mesh_cell_base_3`.

\tparam Gt is the geometric traits class. 
It has to be a model of the concept `RegularTriangulationTraits_3`. 

\tparam MD provides the types of indices used to identify 
the faces of the input complex. It has to be a model 
of the concept `MeshDomain_3`. 

\tparam Tds is the triangulation data structure class to which cells
belong. That parameter is only used by the rebind mechanism (see
`::TriangulationDSCellBase_3::Rebind_TDS`). Users should always use the
default parameter value `void`.

\cgalModels `MeshCellBase_3`

\sa `CGAL::Mesh_complex_3_in_triangulation_3<Tr,CornerIndex,CurveSegmentIndex>` 
\sa `CGAL::Mesh_cell_base_3<Gt, MD, Cb>`

*/
template< typename Gt,  typename MD, typename Tds = void >
class Compact_mesh_cell_base_3 {
public:

/// @}

}; /* end Compact_mesh_cell_base_3 */
} /* end namespace CGAL */
