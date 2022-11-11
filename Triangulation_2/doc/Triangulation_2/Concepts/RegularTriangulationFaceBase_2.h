
/*!
\ingroup PkgTriangulation2Concepts
\cgalConcept

The regular triangulation of a set of weighted points does not
necessarily
have one vertex for each of the input points. Some of the input
weighted points have no cell in the dual power diagrams
and therefore do not correspond to a vertex of the regular
triangulation.
Those weighted points are said to be <I>hidden</I> points.
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

The base face of a regular triangulation
has to be a model
of the concept `RegularTriangulationFaceBase_2` , which refines
the concept `TriangulationFaceBase_2` by adding
in the face a list to store hidden vertices.

\cgalRefines `TriangulationFaceBase_2`

\cgalHasModel `CGAL::Regular_triangulation_face_base_2`

\sa `TriangulationFaceBase_2`
\sa `RegularTriangulationVertexBase_2`

*/

class RegularTriangulationFaceBase_2 {
public:

/// \name Types
/// @{

/*!
A `std::list` of hidden vertices.
*/
typedef std::list<Vertex_handle> Vertex_list;

/// @}

/// \name Access Functions
/// @{

/*!
Returns a reference to the list of vertices hidden by the face.
*/
Vertex_list& vertex_list();

/// @}

}; /* end RegularTriangulationFaceBase_2 */

