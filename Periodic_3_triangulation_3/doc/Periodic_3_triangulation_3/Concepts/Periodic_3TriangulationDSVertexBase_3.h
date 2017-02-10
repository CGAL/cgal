
/*!
\ingroup PkgPeriodic3Triangulation3Concepts
\cgalConcept

A refinement of the concept `TriangulationDSVertexBase_3`
which adds an API for offset.

At the base level of 3D-triangulations
(see Sections \ref P3Triangulation3secdesign and \ref TDS3secdesign),
a vertex provides access to one of its incident cells through a handle.

\cgalRefines `TriangulationDSVertexBase_3`

\cgalHasModel CGAL::Periodic_3_triangulation_ds_vertex_base_3

\sa `TriangulationDataStructure_3`
\sa `TriangulationDSVertexBase_3`
\sa `Periodic_3TriangulationDSCellBase_3`

*/

class Periodic_3TriangulationDSVertexBase_3 {
public:

/// \name Types
/// @{

/*!
A model of the concept `Periodic_3Offset_3`
*/
typedef unspecified_type Periodic_3_offset_3;

/// @}

/// \name Access Functions
/// @{

/*!
Returns the offset stored in the vertex.
*/
Periodic_3_offset_3 offset() const;

/*!
Returns `true` if the offset has been set, `false` otherwise.
*/
bool get_offset_flag() const;

/*!
Sets the offset and sets the offset flag to `true`.
*/
void set_offset(Periodic_3_offset_3 o);

/*!
Sets the offset flag to `false` and clears the offset.
*/
void clear_offset();

/// @}

}; /* end Periodic_3TriangulationDSVertexBase_3 */

