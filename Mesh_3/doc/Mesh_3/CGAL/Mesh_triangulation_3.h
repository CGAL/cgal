namespace CGAL {

/*!
\ingroup PkgMesh_3MeshClasses

The class `Mesh_triangulation_3` is a metafunctor which provides the triangulation to be used as the 3D 
triangulation of a mesh generation process. 

\tparam MD stands for a model of `MeshDomain_3`. 

\tparam Gt stands for a model of `RegularTriangulationTraits_3` 
and defaults to `Kernel_traits<MD>::Kernel`. 

\tparam Concurrency_tag allows to ask for parallel meshing and optimization
algorithms. Possible values are `CGAL::Sequential_tag` (the default) and
`CGAL::Parallel_tag`.

\tparam Vertex_base stands for a model of `MeshVertexBase_3` 
and defaults to `Mesh_vertex_base_3<Gt, MD>`. 

\tparam Cell_base stands for a model of `MeshCellBase_3`
and defaults to `Mesh_cell_base_3<Gt, MD>`. 

\sa `CGAL::make_mesh_3()` 
\sa `CGAL::Mesh_complex_3_in_triangulation_3<Tr,CornerIndex,CurveSegmentIndex>` 

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
`CGAL::Regular_triangulation_3` type 
whose vertex and cell classes are respectively 
`Vertex_base` and `Cell_base`. 
*/ 
typedef Hidden_type type; 

/// @}

}; /* end Mesh_triangulation_3 */
} /* end namespace CGAL */
