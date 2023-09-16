
/*!
\ingroup PkgTriangulation3Concepts
\cgalConcept

The base cell of a Delaunay triangulation must be a model
of the concept `DelaunayTriangulationCellBase_3`, which refines
the concept `TriangulationCellBase_3` by adding
in the cell an operator that computes its circumcenter.


\cgalRefines{TriangulationCellBase_3}

\cgalHasModelsBegin
\cgalHasModels{CGAL::Delaunay_triangulation_cell_base_3}
\cgalHasModels{CGAL::Delaunay_triangulation_cell_base_with_circumcenter_3}
\cgalHasModelsEnd

\sa `DelaunayTriangulationTraits_3`

*/

class DelaunayTriangulationCellBase_3 {
public:

/// \name Types
/// @{
typedef DelaunayTriangulationTraits_3::Point_3 Point;
/// @}

/// \name Access Functions
/// @{
/*!
Returns the circumcenter of the cell.
`DelaunayTriangulationTraits_3` is the geometric traits class of the triangulation.

This operator is required only when the dual functions are called.
*/
const Point& circumcenter(DelaunayTriangulationTraits_3&gt = DelaunayTriangulationTraits_3()) const;
/// @}


}; /* end DelaunayTriangulationCellBase_3 */
