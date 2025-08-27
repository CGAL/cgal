namespace CGAL {

/*!
\ingroup PkgLinearCellComplexClasses

The class `Mesh_data_for_hexmeshing` encapsulates required data for hexmeshing.

\cgalModels{Mesh_data_for_hexmeshing}

\sa `CGAL::Linear_cell_complex_for_combinatorial_map`
\sa `CGAL::Hexmeshing_for_linear_cell_complex_sequential`

*/
struct Mesh_data_for_hexmeshing
{
public:
  /// \name Creation
  /// @{

  /*!
  default constructor
  */
  Mesh_data_for_hexmeshing();
  /*!
  construct with the mesh poly_out
  */
  Mesh_data_for_hexmeshing(Polyhedron_3& poly_out);

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