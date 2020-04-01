
/*!
\ingroup PkgSegmentDelaunayGraph2Concepts
\cgalConcept

The concept `SegmentDelaunayGraphDataStructure_2` refines the
concept `ApolloniusGraphDataStructure_2`. In addition
it provides two methods for the merging of two vertices joined by an
edge of the data structure, and the splitting of a vertex into two.
The method that merges two vertices, called `join_vertices()`
identifies the two vertices and deletes their common two faces. The
method that splits a vertex, called `split_vertex()` introduces a
new vertex that shares an edge and two faces with the old vertex (see
figure below). Notice that the `join_vertices()` and
`split_vertex()` operations are complementary, in the sense that one
reverses the action of the other.

\anchor figsdgdssplitjoin

\cgalFigureBegin{figsdgdssplitjoin,sdg-join_split.png}
The join and split operations. Left to right: The vertex `v` is split
into \f$ v_1\f$ and \f$ v_2\f$. The faces \f$ f\f$ and \f$ g\f$ are
inserted after \f$ f_1\f$ and \f$ f_2\f$, respectively, in the
counter-clockwise sense. The vertices \f$ v_1\f$, \f$ v_2\f$ and the
faces \f$ f\f$ and \f$ g\f$ are returned as a `boost::tuple` in that
order. Right to left: The edge `(f,i)` is collapsed, and thus the
vertices \f$ v_1\f$ and \f$ v_2\f$ are joined. The vertex `v` is
returned.
\cgalFigureEnd

We only describe the additional requirements with respect to the
`ApolloniusGraphDataStructure_2` concept.

\cgalRefines `ApolloniusGraphDataStructure_2`

\cgalHasModel `CGAL::Triangulation_data_structure_2<Vb,Fb>`

\sa `TriangulationDataStructure_2`
\sa `ApolloniusGraphDataStructure_2`
\sa `SegmentDelaunayGraphVertexBase_2`
\sa `TriangulationFaceBase_2`

*/

class SegmentDelaunayGraphDataStructure_2 {
public:

/// \name Modification
/// @{

/*!
Joins
the vertices that are endpoints of the edge `(f,i)` and returns
a vertex handle to common vertex.
*/
Vertex_handle join_vertices(Face_handle f, int i);

/*!
Splits the vertex `v` into two vertices `v1` and
`v2`. The common faces `f` and `g` of `v1` and
`v2` are created after (in the counter-clockwise sense) the
faces `f1` and `f2`. The 4-tuple `(v1,v2,f,g)` is
returned (see \cgalFigureRef{figsdgdssplitjoin}).
*/
boost::tuples::tuple<Vertex_handle, Vertex_handle, Face_handle,
Face_handle>
split_vertex(Vertex_handle v, Face_handle f1, Face_handle
f2);

/// @}

}; /* end SegmentDelaunayGraphDataStructure_2 */

