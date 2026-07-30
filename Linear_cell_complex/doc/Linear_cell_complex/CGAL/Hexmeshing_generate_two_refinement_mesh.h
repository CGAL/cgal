namespace CGAL {

/*!
\ingroup PkgLinearCellComplexConstructions

generates a pure hexahedral mesh from a triangle mesh using the two refinement algorithm described in \cgalCite{cgal:owen2017template-based}.

\tparam TriangleMesh a model of `FaceListGraph`
\param tmesh a triangle mesh 
\param cube_cells_per_dim Grid cells per dimension
\param nb_levels How many times to perform refinement
\param trim Whether to apply trimming to remove excess volumes after refinement. A volume cell is considered to be excess if the full volume of it is outside of the surface mesh.
              (default: `false`)

\sa `CGAL::Linear_cell_complex_for_combinatorial_map`

*/
  template <typename TriangleMesh>
  LCC generate_hexahedral_mesh_using_two_refinement
  (const TriangleMesh& tmesh, int cube_cells_per_dim, int nb_levels, bool trim=false);

} // end namespace CGAL
