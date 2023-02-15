namespace CGAL {

/*!
\ingroup BGLGraphExternalIndices

The class `Triangulation_vertex_base_with_id_2` is a model of the
concept `TriangulationVertexBase_2`, the base vertex of a
2D-triangulation.  It provides an integer field that can be used to
index vertices for \bgl algorithms.

Note that the user is in charge of setting indices correctly before
running a graph algorithm, by calling the function
`CGAL::set_triangulation_ids(Triangulation&)`.

\tparam TriangulationTraits_2 is the geometric traits class
and must be a model of `TriangulationTraits_2`.

\tparam TriangulationVertexBase_2 must be a vertex base class from which
`Triangulation_vertex_base_with_id_2` derives. It has the default
value `Triangulation_vertex_base_2<TriangulationTraits_2>`.

\cgalModels `TriangulationVertexBase_2`

\sa `CGAL::Triangulation_vertex_base_2`
*/
template< typename TriangulationTraits_2, typename TriangulationVertexBase_2 >
class Triangulation_vertex_base_with_id_2 : public TriangulationVertexBase_2 {
public:

/// \name Access Functions
/// @{

/*!
Returns the index.
*/
int id() const;

/*!
Returns a reference to the index stored in the vertex.
*/
int& id();

/// @}

}; /* end Triangulation_vertex_base_with_id_2 */
} /* end namespace CGAL */
