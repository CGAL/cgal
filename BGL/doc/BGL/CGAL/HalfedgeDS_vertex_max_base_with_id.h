
namespace CGAL {

/*!
\ingroup BGLGraphExternalIndices

The class `HalfedgeDS_vertex_max_base_with_id` is a model of the `HalfedgeDSVertex`
concept.  It is
equivalent to `HalfedgeDS_vertex_base< Refs, Tag_true>`
with an added integer field which can be used to index vertices
in \bgl algorithms..
The class contains support for the point and the required type definitions.
It can be used for deriving own vertices.

Note that the user is in charge to set the index correctly before
running a graph algorithm.

\tparam Refs must be an instantiation of a `HalfedgeDS`.

\cgalModels `HalfedgeDSVertex`

\sa `HalfedgeDS<Traits,Items,Alloc>`
\sa `HalfedgeDSItems`
\sa `PolyhedronItems_3`
\sa `CGAL::HalfedgeDS_min_items`
\sa `CGAL::HalfedgeDS_halfedge_min_base<Refs>`
\sa `CGAL::HalfedgeDS_face_min_base<Refs>`
\sa `CGAL::HalfedgeDS_vertex_base<Refs>`
\sa `CGAL::HalfedgeDS_halfedge_max_base_with_id<Refs>`
\sa `CGAL::HalfedgeDS_face_max_base_with_id<Refs>`
\sa `CGAL::Polyhedron_items_with_id_3`

*/
template< typename Refs >
class HalfedgeDS_vertex_max_base_with_id {
public:

/// \name Creation
/// @{

/*!
%Default constructor.
*/
HalfedgeDS_vertex_max_base_with_id();

/*!
Returns the index.
*/
int id() const;

/*!
Returns a reference to the index stored in the vertex.
*/
int& id();

/// @}

}; /* end HalfedgeDS_vertex_max_base_with_id */
} /* end namespace CGAL */
