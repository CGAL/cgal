namespace CGAL {

/*!
\ingroup PkgLinearCellComplexConstructions

Two refinement algorithm described by \cgalCite{cgal:owen2017template-based} on a 3-dimensional linear cell complex.
It performs grid-based hexahedral meshing for a given surface mesh with multiple refinement levels.

\param file Path to an .off surface mesh file
\param cube_cells_per_dim Grid cells per dimension
\param nb_levels How many times to perform refinement
\param trim Whether to apply trimming to remove excess volumes after refinement. A volume cell is considered to be excess if the full volume of it is outside of the surface mesh.
              (default: false)

\sa `CGAL::Linear_cell_complex_for_combinatorial_map`

*/
  LCC generate_two_refinement_mesh(const std::string& file, int cube_cells_per_dim, int nb_levels, bool trim=false);

/*!
\ingroup PkgLinearCellComplexConstructions

Two refinement algorithm described by \cgalCite{cgal:owen2017template-based} on a 3-dimensional linear cell complex.
It performs grid-based hexahedral meshing for a given surface mesh with multiple refinement levels.

\tparam TriangleMesh model of FaceGraph
              (default: Polyhedron)
\param poly Surface mesh
\param cube_cells_per_dim Grid cells per dimension
\param nb_levels How many times to perform refinement
\param trim Whether to apply trimming to remove excess volumes after refinement. A volume cell is considered to be excess if the full volume of it is outside of the surface mesh.
              (default: false)

\sa `CGAL::Linear_cell_complex_for_combinatorial_map`

*/
  template <typename TriangleMesh=Polyhedron>
  LCC generate_two_refinement_mesh(TriangleMesh& poly, int cube_cells_per_dim, int nb_levels, bool trim=false);

} // end namespace CGAL