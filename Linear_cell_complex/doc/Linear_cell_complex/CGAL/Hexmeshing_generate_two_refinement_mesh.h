namespace CGAL {

/*!
\ingroup PkgLinearCellComplexConstructions

generates a pure hexahedral mesh from a triangle mesh using the two-refinement algorithm described in \cgalCite{cgal:owen2017template-based}.
Starts to create a regular grid of `cube_cells_per_dim`\f$ ^3\f$ voxels. Then refines voxels intersected by the surface `nb_levels` times, while creating transitions between refined and non-refined voxels.

\tparam TriangleMesh a model of `FaceListGraph`
\param tmesh a triangle mesh
\param cube_cells_per_dim number of grid cells, per dimension
\param nb_levels How many times to perform refinement
\param trim `true` to apply trimming, i.e., remove volumes that are entirely outside of the surface mesh
\param smooth `true` to smooth the hexahedral mesh, using a 3D Laplacian smoothing method

\sa `CGAL::Linear_cell_complex_for_combinatorial_map`

*/
  template <typename TriangleMesh>
  LCC generate_hexahedral_mesh_using_two_refinement
  (const TriangleMesh& tmesh, int cube_cells_per_dim, int nb_levels, bool trim=true, bool smooth=true);

} // end namespace CGAL
