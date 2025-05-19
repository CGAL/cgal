/*!
\ingroup PkgConstrainedTriangulation3Concepts
\cgalConcept

The concept `ConformingConstrainedDelaunayTriangulationCellBase_3` refines the concept
`TriangulationCellBase_3` and and describes the requirements for a base cell class of
the `CGAL::Conforming_constrained_Delaunay_triangulation_3` class.

\cgalRefines{TriangulationCellBase_3, BaseWithTimeStamp}

\cgalHasModelsBegin
\cgalHasModels{CGAL::Conforming_constrained_Delaunay_triangulation_cell_base_3}
\cgalHasModelsEnd

\sa `ConformingConstrainedDelaunayTriangulationVertexBase_3`
*/
class ConformingConstrainedDelaunayTriangulationCellBase_3 {
public:

  /// @name Access Functions
  ///
  /// The following functions return a reference to an object of type
  /// `CGAL::Conforming_constrained_Delaunay_triangulation_cell_data_3`, that contains
  /// the per-cell data required by the implementation of the
  /// `CGAL::Conforming_constrained_Delaunay_triangulation_3` class.
  /// @{
  CGAL::Conforming_constrained_Delaunay_triangulation_cell_data_3& ccdt_3_data();
  const CGAL::Conforming_constrained_Delaunay_triangulation_cell_data_3& ccdt_3_data() const;
  /// @}
};
