namespace CGAL {

/*!
\ingroup PkgHexmeshingLinearCellComplex

The class `MeshDataForHexmeshing` encapsulates required data for hexmeshing.

\cgalModels{MeshDataForHexmeshing}

\sa `CGAL::Linear_cell_complex_for_combinatorial_map`

*/
struct MeshDataForHexmeshing
{
public:
  /// \name Creation
  /// @{

  /*!
  default constructor
  */
  MeshDataForHexmeshing();

  /// @}

  /// \name Operations
  /// @{

  /*!
  imports polygon mesh data from .off file.
  */
  void load_surface(const std::string& file);

  /*!
  construct cube_cells_per_dim^3 grid from the imported polygon mesh.
  */
  void cubic_grid_from_aabb(int cube_cells_per_dim);

  /// @}
};

} // end namespace CGAL