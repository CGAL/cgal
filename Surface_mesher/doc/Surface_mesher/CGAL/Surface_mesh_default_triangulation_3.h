
namespace CGAL {

/*!
\ingroup PkgSurfaceMesher3Classes

The class `Surface_mesh_default_triangulation_3` is a model of the concept 
`SurfaceMeshTriangulation_3`, whose vertex base and cell base classes 
are models of the concepts 
`SurfaceMeshVertexBase_3` and 
`SurfaceMeshCellBase_3` respectively. 

\cgalModels `SurfaceMeshTriangulation_3`

\sa `Surface_mesh_complex_2_in_triangulation_3<Tr>` 
\sa `make_surface_mesh` 

*/

class Surface_mesh_default_triangulation_3 {
public:
  
  /// the vertex type which is model of `::SurfaceMeshVertexBase_3`.
  typedef unspecified_type Vertex;


  /// the cell type which is model of `::SurfaceMeshCellBase_3`.
  typedef unspecified_type Cell;  

/// @}

}; /* end Surface_mesh_default_triangulation_3 */
} /* end namespace CGAL */
