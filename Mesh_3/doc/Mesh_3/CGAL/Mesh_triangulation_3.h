namespace CGAL {

/*!
\ingroup PkgMesh_3MeshClasses

The class `Mesh_triangulation_3` provides a default triangulation to be used as the 3D 
triangulation of a mesh generation process. 

Parameters 
-------------- 

Template parameter `MD` stands for a model of `MeshDomain_3`. 

Template parameter `Gt` stands for a model of `RegularTriangulationTraits_3` 
and defaults to `Kernel_traits<MD>::Kernel`. 

\sa `CGAL::make_mesh_3` 
\sa `CGAL::Mesh_complex_3_in_triangulation_3<Tr,CornerIndex,CurveSegmentIndex>` 

*/
template< typename MD, typename Gt >
class Mesh_triangulation_3 {
public:

/// \name Types 
/// @{

/*! 
`CGAL::Regular_triangulation_3` type 
whose vertex and cell classes are respectively 
`Mesh_vertex_base_3<MD,Gt>` and 
`Mesh_cell_base_3<MD,Gt>`. 
*/ 
typedef Hidden_type type; 

/// @}

}; /* end Mesh_triangulation_3 */
} /* end namespace CGAL */
