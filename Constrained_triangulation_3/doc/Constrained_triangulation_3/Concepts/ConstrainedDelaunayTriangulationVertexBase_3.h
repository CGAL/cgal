/*!
\cgalConcept

The concept `ConstrainedDelaunayTriangulationVertexBase_3` refines the concept
`TriangulationVertexBase_3` and is the base vertex class for
the `CGAL::Constrained_Delaunay_triangulation_3` class template.

\cgalRefines{TriangulationVertexBase_3}

\cgalHasModelsBegin
\cgalHasModels{CGAL::Constrained_Delaunay_triangulation_vertex_base_3}
\cgalHasModelsEnd


\sa `ConstrainedDelaunayTriangulationCellBase_3`

\todo Verify that the See also section is correctly formatted by Doxygen.
*/
class ConstrainedDelaunayTriangulationVertexBase_3 {
public:

  /// @name Access Functions
  ///
  /// The following functions provide a reference to an object of type
  /// `CGAL::Constrained_Delaunay_triangulation_vertex_data_3`, that contains
  /// the per-vertex data required by the implementation of the
  /// `CGAL::Constrained_Delaunay_triangulation_3` class template.
  /// @{
  CGAL::Constrained_Delaunay_triangulation_vertex_data_3& cdt_3_data();
  const CGAL::Constrained_Delaunay_triangulation_vertex_data_3& cdt_3_data() const;
  /// @}
}
