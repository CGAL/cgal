
/*!
\ingroup PkgTriangulation2Concepts
\cgalConcept

The vertex of a triangulation
included in a triangulation hierarchy has to provide
some pointers to the corresponding vertices in the
triangulations of the next and preceding levels.
Therefore, the concept `TriangulationHierarchyVertexBase_2`
refines the concept `TriangulationVertexBase_2`,
adding handles to the corresponding vertices in the
next and previous level triangulations.

\cgalRefines{TriangulationVertexBase_2}

\cgalHasModelsBegin
\cgalHasModels{CGAL::Triangulation_hierarchy_vertex_base_2<Vb>}
\cgalHasModelsEnd

\sa `CGAL::Triangulation_hierarchy_2<Tr>`

*/

class TriangulationHierarchyVertexBase_2 {
public:

/// \name Operations
/// @{

/*!
returns the corresponding vertex
(if any) of the next level triangulation.
*/
Vertex_handle up();

/*!
returns the corresponding vertex
of the previous level triangulation.
*/
Vertex_handle down();

/*!
sets the handle pointing to to the corresponding vertex
of the next level triangulation.
*/
void set_up(Vertex_handle u);

/*!
sets the handle pointing to the corresponding vertex
of the previous level triangulation.
*/
void set_down(Vertex_handle d);

/// @}

}; /* end TriangulationHierarchyVertexBase_2 */

