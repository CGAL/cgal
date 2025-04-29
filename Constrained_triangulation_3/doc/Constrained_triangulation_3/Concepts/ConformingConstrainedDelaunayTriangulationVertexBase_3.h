/*!
\ingroup PkgConstrainedTriangulation3Concepts
\cgalConcept

The concept `ConformingConstrainedDelaunayTriangulationVertexBase_3` refines the concept
`TriangulationVertexBase_3` and is the base vertex class for
the triangulation returned by
`CGAL::make_conforming_constrained_Delaunay_triangulation_3()` function template.

\cgalRefines{TriangulationVertexBase_3, BaseWithTimeStamp}

\cgalHasModelsBegin
\cgalHasModels{CGAL::Conforming_constrained_Delaunay_triangulation_vertex_base_3}
\cgalHasModelsEnd

\sa `ConformingConstrainedDelaunayTriangulationCellBase_3`

*/
class ConformingConstrainedDelaunayTriangulationVertexBase_3 {
public:

  /// @name Access Functions
  ///
  /// The following functions return a reference to an object of type
  /// `CGAL::Conforming_constrained_Delaunay_triangulation_vertex_data_3`, that contains
  /// the per-vertex data required by the implementation of the
  /// `CGAL::make_conforming_constrained_Delaunay_triangulation_3()` function template.
  /// @{
  CGAL::Conforming_constrained_Delaunay_triangulation_vertex_data_3& ccdt_3_data();
  const CGAL::Conforming_constrained_Delaunay_triangulation_vertex_data_3& ccdt_3_data() const;
  /// @}
}
