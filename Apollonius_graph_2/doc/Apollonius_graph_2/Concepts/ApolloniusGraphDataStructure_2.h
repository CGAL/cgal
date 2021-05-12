
/*!
\ingroup PkgApolloniusGraph2Concepts
\cgalConcept

The concept `ApolloniusGraphDataStructure_2` refines the concept
`TriangulationDataStructure_2`. In addition
it provides two methods for the insertion and removal of a degree 2
vertex in the data structure. The insertion method adds a new vertex
to the specified edge, thus creating two new edges. Moreover, it
creates two new faces that have the two newly created edges in
common (see figure below). The removal method performs the reverse
operation.

\image html insert_degree_2.png
\image latex insert_degree_2.png

<center><b>Insertion and removal of degree 2 vertices. Left to right:
The edge `(f,i)` is replaced by two edges by means of inserting a
vertex `v` on the edge. The faces \f$ f_1\f$ and \f$ f_2\f$ are
created. Right to left: the faces \f$ f_1\f$ and \f$ f_2\f$ are
destroyed. The vertex `v` is deleted and its two adjacent edges are
merged. </b></center>

We only describe the additional requirements with respect to the
`TriangulationDataStructure_2` concept.

\cgalRefines `TriangulationDataStructure_2`

\cgalHasModel `CGAL::Triangulation_data_structure_2<Vb,Fb>`

\sa `TriangulationDataStructure_2`
\sa `ApolloniusGraphVertexBase_2`
\sa `TriangulationFaceBase_2`

*/

class ApolloniusGraphDataStructure_2 {
public:

/// \name Insertion
/// @{

/*!
inserts a degree two vertex and two faces adjacent to it that have two common edges.

The edge defined by the face handle `f` and the integer `i` is duplicated. It returns a handle
to the vertex created.
*/
Vertex_handle insert_degree_2(Face_handle f, int i);

/// @}

/// \name Removal
/// @{

/*!
Removes a degree 2
vertex and the two faces adjacent to it. The two edges of the star of
`v` that are not incident to it are collapsed.
\pre The degree of `v` must be equal to 2.
*/
void remove_degree_2(Vertex_handle v);

/// @}

}; /* end ApolloniusGraphDataStructure_2 */

