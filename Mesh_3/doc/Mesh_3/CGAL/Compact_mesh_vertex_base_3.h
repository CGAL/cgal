namespace CGAL {

/*!
\ingroup PkgMesh_3MeshClasses

The class `Compact_mesh_vertex_base_3` is a model of the concept `MeshVertexBase_3`. 
It is designed to serve as vertex base class for the 3D triangulation 
used in a 3D mesh generation process. It is more compact in memory than
`Mesh_vertex_base_3`.

\tparam Gt is the geometric traits class. 
It must be a model of the concept `RegularTriangulationTraits_3`. 

\tparam MD provides the types of indices 
used to identify 
the faces of the input complex. It must be a model 
of the concept `MeshDomain_3`. 
	
\tparam Vb is the vertex base class. It has to be a model 
of the concept `TriangulationVertexBase_3` and defaults to 
`Triangulation_vertex_base_3<Gt>`. 

\cgalModels `MeshVertexBase_3` 

\sa `CGAL::Mesh_complex_3_in_triangulation_3<Tr,CornerIndex,CurveSegmentIndex>` 

*/
template< typename Gt, typename MD, typename Vb >
class Compact_mesh_vertex_base_3 : public Vb {
public:

/// @}

}; /* end Compact_mesh_vertex_base_3 */
} /* end namespace CGAL */
