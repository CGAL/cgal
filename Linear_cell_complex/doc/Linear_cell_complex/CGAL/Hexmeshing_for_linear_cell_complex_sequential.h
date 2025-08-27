namespace CGAL {

/*!
\ingroup PkgLinearCellComplexClasses

The class `Hexmeshing_for_linear_cell_complex` contains two refinement algorithm described by \cgalCite{cgal:owen2017template-based} on a 3-dimensional linear cell complex.

\cgalModels{Hexmeshing_for_linear_cell_complex}

\sa `CGAL::Linear_cell_complex_for_combinatorial_map`
\sa `CGAL::Hexmeshing_mesh_data_for_hexmeshing`

*/
struct Hexmeshing_for_linear_cell_complex
{
public:
  /// \name Creation
  /// @{

  /*!
  default constructor
  */
  Hexmeshing_for_linear_cell_complex();

  /// @}

  /// \name Operations
  /// @{

  /*!
  This function implements the two_refinement algorithm as described in the \cgalCite{cgal:owen2017template-based},
  operating on a 3-dimensional linear cell complex. It performs grid-based hexahedral meshing
  with multiple refinement levels.
  \pre mesh have polygon mesh data as well as grid data.
  */
  void two_refinement_algorithm(Mesh_data_for_hexmeshing &mesh, int nb_levels, bool trim);

  /// @}
};

} // end namespace CGAL