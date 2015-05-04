#if 0
namespace boost {

/*!
\ingroup PkgBGLTraits  

The class `CGAL::Surface_mesh` is an undirected graph, does not
allow multiedges and models the following concepts:

\anchor GraphTraitsSurfaceMesh

- \bgllink{BidirectionalGraph}
- \bgllink{`VertexAndEdgeListGraph`}
- \bgllink{PropertyGraph}
- `HalfedgeGraph`
- `HalfedgeListGraph`
- `MutableHalfedgeGraph`
- `FaceGraph`
- `FaceListGraph`
- `MutableFaceGraph`

All `.._size_type` are equal to the `size_type` of `CGAL::Surface_mesh`.

The following interior properties are always supported:

Enum                   |  PropertyMap model
---------------------- | -------------------
`vertex_index_t`       |  ReadablePropertyMap
`halfedge_index_t`     |  ReadablePropertyMap
`edge_index_t`         |  ReadablePropertyMap
`face_index_t`         |  ReadablePropertyMap
`halfedge_is_border_t` |  ReadablePropertyMap
`vertex_is_border_t`   |  ReadablePropertyMap
`vertex_point_t`       |  LvaluePropertyMap


`CGAL::Surface_mesh` comes with its own dynamic property mechanism to
extend it with other properties on the fly.

\attention `CGAL::Surface_mesh` does not model \bgllink{MutableGraph}
as operations such as `add_edge(v,w,g)` cannot be implemented. 
None of the mutating operations preserves invariants of the halfedge
data structure and care has to be taken to maintain validity when
those operations are used.
*/
template<>
class graph_traits< CGAL::Surface_mesh > {}; 
} /* end namespace boost */
#endif
