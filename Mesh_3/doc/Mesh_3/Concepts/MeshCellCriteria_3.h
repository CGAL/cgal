/*!
\ingroup PkgMesh3Concepts
\cgalConcept

The Delaunay refinement process involved in the
template functions `CGAL::make_mesh_3()` and `CGAL::refine_mesh_3()`
is guided by a set of elementary refinement criteria
that concern either mesh tetrahedra or surface facets.
The concept `MeshCellCriteria_3` describes criteria for mesh tetrahedra.
For a mesh complex type `C3T3`, the criteria are evaluated on cells of its
nested triangulation.

\cgalHasModelsBegin
\cgalHasModels{CGAL::Mesh_cell_criteria_3<C3T3>}
\cgalHasModelsEnd

\sa `MeshEdgeCriteria_3`
\sa `MeshFacetCriteria_3`
\sa `MeshCriteria_3`
\sa `CGAL::make_mesh_3()`
\sa `CGAL::refine_mesh_3()`

*/

class MeshCellCriteria_3 {
public:

/// \name Types
/// @{

/*!
Type representing the quality of a
cell. Must be a model of CopyConstructible and
LessThanComparable. Between two cells, the one which has the lower
quality must have the lower `Cell_quality`.
*/
typedef unspecified_type Cell_quality;

/*!
Type representing if a cell is bad or not. Must
be contextually convertible to `bool`. If it converts to `true` then the cell
is bad, otherwise the cell is good with regard to the criteria.

In addition, an object of this type must contain an object of type
`Cell_quality` if it represents
a bad cell. `Cell_quality` must be accessible by `operator*()`.
Note that `std::optional<Cell_quality>` is a natural model of this concept.
*/
typedef unspecified_type Is_cell_bad;

/// @}

/// \name Operations
/// @{

/*!
Returns the `Is_cell_bad` value for the cell `c` of the triangulation nested
in `c3t3`, the mesh complex used by the mesh generation function.
*/
Is_cell_bad operator()(const C3T3& c3t3,
                       typename C3T3::Triangulation::Cell_handle c) const;

/// @}

}; /* end MeshCellCriteria_3 */
