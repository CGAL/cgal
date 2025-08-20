namespace CGAL {

/*!
\ingroup PkgHexmeshingLinearCellComplex

The class `HexMeshingData` contains two refinement algorithm described by [2017 論文] on a 3-dimensional linear cell complex.

\cgalModels{HexMeshingData}

\sa `CGAL::Linear_cell_complex_for_combinatorial_map`

*/
struct HexMeshingData
{
public:
  /// \name Creation
  /// @{

  /*!
  default constructor
  */
  HexMeshingData();

  /// @}

  /// \name Operations
  /// @{

  /*!
  This function implements the two_refinement algorithm as described in the [2017 paper],
  operating on a 3-dimensional linear cell complex. It performs grid-based hexahedral meshing
  with multiple refinement levels.
  \pre mesh have polygon mesh data as well as grid data.
  */
  void two_refinement_algorithm(MeshDataForHexmeshing &mesh, int nb_levels, bool trim);

  /// @}
};

} // end namespace CGAL