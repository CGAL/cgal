/*!
\ingroup PkgCT_3Concepts
\cgalConcept

The concept `ConstrainedDelaunayTriangulationCellBase_3` refines the concept
`TriangulationCellBase_3` and is the base cell class for
the `CGAL::make_constrained_Delaunay_triangulation_3()` function template.

\cgalRefines{TriangulationCellBase_3}

\cgalHasModelsBegin
\cgalHasModels{CGAL::Constrained_Delaunay_triangulation_cell_base_3}
\cgalHasModelsEnd

\sa `ConstrainedDelaunayTriangulationVertexBase_3`
*/
class ConstrainedDelaunayTriangulationCellBase_3 {
public:

  /// @name Access Functions
  ///
  /// The following functions return a reference to an object of type
  /// `CGAL::Constrained_Delaunay_triangulation_cell_data_3`, that contains
  /// the per-cell data required by the implementation of the
  /// `CGAL::make_constrained_Delaunay_triangulation_3()` function template.
  /// @{
  CGAL::Constrained_Delaunay_triangulation_cell_data_3& cdt_3_data();
  const CGAL::Constrained_Delaunay_triangulation_cell_data_3& cdt_3_data() const;
  /// @}
};
