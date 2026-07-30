namespace CGAL {

/*!
\ingroup PkgLinearCellComplexConstructions

Two refinement algorithm described by \cgalCite{cgal:owen2017template-based} on a 3-dimensional linear cell complex.
It performs grid-based hexahedral meshing for a given surface mesh with multiple refinement levels.

\tparam TriangleMesh a model of FaceListGraph
\param tmesh a triangle mesh 
\param cube_cells_per_dim Grid cells per dimension
\param nb_levels How many times to perform refinement
\param trim Whether to apply trimming to remove excess volumes after refinement. A volume cell is considered to be excess if the full volume of it is outside of the surface mesh.
              (default: false)

\sa `CGAL::Linear_cell_complex_for_combinatorial_map`

*/
  template <typename TriangleMesh>
  LCC generate_hexahedral_mesh_using_two_refinement
  (const TriangleMesh& tmesh, int cube_cells_per_dim, int nb_levels, bool trim=false);

} // end namespace CGAL
