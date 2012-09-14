namespace CGAL {

/*!
\ingroup PkgMesh_3MeshClasses

The class `Mesh_cell_base_3` is a model of the concept `MeshCellBase_3`. 
It is designed to serve as cell base class for the 3D triangulation 
used in the 3D mesh generation process. 

\tparam MD provides the types of indices used to identify 
the faces of the input complex. It has to be a model 
of the concept `MeshDomain_3`. 

\tparam Gt is the geometric traits class. 
It has to be a model of the concept `RegularTriangulationTraits_3`. 

\tparam Cb is the cell base class. It has to be a model 
of the concept `RegularTriangulationCellBase_3` and defaults to 
`Regular_triangulation_cell_base_3<Gt>`. 

\models ::MeshCellBase_3 

\sa `MeshCellBase_3` 
\sa `CGAL::Mesh_complex_3_in_triangulation_3<Tr,CornerIndex,CurveSegmentIndex>` 

*/
template< typename MD, typename Gt, typename Cb >
class Mesh_cell_base_3 : public Cb {
public:

/// @}

}; /* end Mesh_cell_base_3 */
} /* end namespace CGAL */
