
/*!
\ingroup PkgPeriodic3Triangulation3Concepts
\cgalConcept

A refinement of the concept `TriangulationDSCellBase_3`
which adds an API for offsets.

At the base level (see Sections \ref P3Triangulation3secdesign
and \ref TDS3secdesign), a cell stores handles to its four vertices
and to its four neighbor cells. The vertices and neighbors are
indexed 0, 1, 2 and 3. Neighbor `i` lies opposite to vertex `i`.

For periodic triangulations, the cell base class needs to
additionally store an offset for each vertex. Only the last three
bits of each integer are required to be stored. The remaining part
does not contain any information.

\cgalRefines `TriangulationDSCellBase_3`

\cgalHasModel CGAL::Periodic_3_triangulation_ds_cell_base_3

\sa `TriangulationDataStructure_3`
\sa `TriangulationDSCellBase_3`
\sa `Periodic_3TriangulationDSVertexBase_3`

*/

class Periodic_3TriangulationDSCellBase_3 {
public:

/// \name Access Functions
/// @{

/*!
Returns the offset of vertex `i`.
\pre \f$ i \in\{0, 1, 2, 3\}\f$.
*/
int offset(int i) const;

/// @}

/// \name Setting
/// @{

/*!
Sets the vertex offsets according to `off0` to `off3`.
*/
void set_offsets(int off0, int off1, int off2, int off3);

/// @}

}; /* end Periodic_3TriangulationDSCellBase_3 */

