
/*!
\ingroup PkgTriangulation2Concepts
\cgalConcept

The regular triangulation of a set of weighted points does not
necessarily
have one vertex for each of the input points. Some of the input
weighted points have no cell in the dual power diagrams
and therefore do not correspond to a vertex of the regular
triangulation.
Those weighted point are said to be <I>hidden</I> points.
A point which is hidden at a given time may appear later as a vertex of
the regular triangulation upon removal on some other weighted point.
Therefore, hidden points have to be stored somewhere.
The regular triangulation store those hidden points
in special vertices called <I>hidden</I> vertices.

A hidden point can appear as vertex of the triangulation
only when the
two dimensional face where its point component is located
(the face which hides it)
is removed. Therefore we decided to store
in each face of a regular triangulation
the list of hidden vertices whose points are located in the face.
Thus points hidden by a face are easily reinserted in the triangulation
when the face is removed.

The base vertex of a regular triangulation has to be a model
of the concept `RegularTriangulationVertexBase_2` .
The concept `RegularTriangulationVertexBase_2` refines the concept
`TriangulationVertexBase_2`,
just adding a Boolean to mark if the vertex is a
vertex of the triangulation or a hidden vertex.

\cgalRefines `TriangulationVertexBase_2`

\cgalHasModel `CGAL::Regular_triangulation_vertex_base_2`

\sa `TriangulationVertexBase_2`
\sa `CGAL::Regular_triangulation_vertex_base_2`

*/

class RegularTriangulationVertexBase_2 {
public:

/*!
Must be the same as the point type `RegularTriangulationTraits_2::Weighted_point_2`
defined by the geometric traits class of the triangulation.
*/
typedef unspecified_type Point;

/// \name Access Functions
/// @{

/*!
returns `true`, iff the vertex is hidden.
*/
bool is_hidden();

/*!
Mark the vertex as hidden or as not hidden.
*/
void set_hidden(bool b);

/// @}

}; /* end RegularTriangulationVertexBase_2 */

