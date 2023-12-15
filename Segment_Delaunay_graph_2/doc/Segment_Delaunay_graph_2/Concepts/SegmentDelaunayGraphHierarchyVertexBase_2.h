
/*!
\ingroup PkgSegmentDelaunayGraph2Concepts
\cgalConcept

The vertex of a segment Delaunay graph
included in a segment Delaunay graph hierarchy has to provide
some pointers to the corresponding vertices in the
graphs of the next and preceding levels.
Therefore, the concept `SegmentDelaunayGraphHierarchyVertexBase_2`
refines the concept `SegmentDelaunayGraphVertexBase_2`, by
adding two vertex handles to the corresponding vertices for the
next and previous level graphs.

\cgalRefines{SegmentDelaunayGraphVertexBase_2}

\cgalHeading{Types}

`SegmentDelaunayGraphHierarchyVertexBase_2` does not introduce
any types in addition to those of
`SegmentDelaunayGraphVertexBase_2`.

\cgalHeading{Creation}

The `SegmentDelaunayGraphHierarchyVertexBase_2` concept does not
introduce any constructors in addition to those of the
`SegmentDelaunayGraphVertexBase_2` concept.

\cgalHasModelsBegin
\cgalHasModels{CGAL::Segment_Delaunay_graph_hierarchy_vertex_base_2<CGAL::Segment_Delaunay_graph_vertex_base_2<St> >}
\cgalHasModelsEnd

\sa `SegmentDelaunayGraphDataStructure_2`
\sa `SegmentDelaunayGraphVertexBase_2`
\sa `CGAL::Segment_Delaunay_graph_hierarchy_2<Gt,St,SSTag,DS>`
\sa `CGAL::Triangulation_data_structure_2<Vb,Fb>`
\sa `CGAL::Segment_Delaunay_graph_vertex_base_2<St,Vb>`
\sa `CGAL::Segment_Delaunay_graph_hierarchy_vertex_base_2<Vbb>`

*/

class SegmentDelaunayGraphHierarchyVertexBase_2 {
public:

/// \name Operations
/// @{

/*!
Returns a handle to the corresponding
vertex of the next level segment Delaunay graph. If such a vertex
does not exist `Vertex_handle()` is returned.
*/
Vertex_handle up();

/*!
Returns a handle to the corresponding
vertex of the previous level segment Delaunay graph. If such a
vertex does not exist `Vertex_handle()` is returned.
*/
Vertex_handle down();

/*!
Sets the handle for the
vertex of the next level segment Delaunay graph.
*/
void set_up(Vertex_handle u);

/*!
Sets the handle for the
vertex of the previous level segment Delaunay graph.
*/
void set_down(Vertex_handle d);

/// @}

}; /* end SegmentDelaunayGraphHierarchyVertexBase_2 */

