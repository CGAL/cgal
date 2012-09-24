namespace CGAL {

/*!
\ingroup PkgBGLEnums
The constant `edge_is_border` is a 
<A HREF="http://www.boost.org/libs/graph/doc/PropertyTag.html">property tag</A> 
which refers to the property of an edge of being a border edge.

`edge_is_border` is an 
<A HREF="http://www.boost.org/libs/graph/doc/using_property_maps.html">interior property</A>. 
That is, a
<A HREF="http://www.boost.org/libs/property_map/doc/property_map.html">property map</A> 
for `edge_is_border` can be extracted from any model of a `HalfedgeGraph`
using the <span class="textsc">Bgl</span>
<A HREF="http://www.boost.org/libs/graph/doc/PropertyGraph.html">PropertyGraph</A> interface: `boost::get(edge_is_border,graph)`

The Boolean flag that indicates if the edge is a border can be directly accessed via:
`boost::get(edge_is_border,graph,edge)`.
*/
enum edge_is_border_t { edge_is_border };
/*!
\ingroup PkgBGLEnums

The constant `vertex_is_border` is a 
<A HREF="http://www.boost.org/libs/graph/doc/PropertyTag.html">property tag</A> which refers to the property
of a vertex of being a border vertex.

`vertex_is_border` is an 
<A HREF="http://www.boost.org/libs/graph/doc/using_property_maps.html">interior property</A>, that is, a
<A HREF="http://www.boost.org/libs/property_map/doc/property_map.html">property map</A> 
for `vertex_is_border` can be extracted from any  model of 
a `HalfedgeGraph` using the <span class="textsc">Bgl</span>
<A HREF="http://www.boost.org/libs/graph/doc/PropertyGraph.html">PropertyGraph</A> interface:
`boost::get(vertex_is_border,graph)`

The Boolean flag that indicates if the vertex is a border can be directly accessed via:
`boost::get(vertex_is_border,graph,edge)`
*/
enum vertex_is_border_t { vertex_is_border };

/*!
\ingroup PkgBGLEnums
The constant `vertex_point` is a <A HREF="http://www.boost.org/libs/graph/doc/PropertyTag.html">property tag</A> which refers to the  geometric embedding property of a  vertex of a `HalfedgeGraph`.

A `vertex_point` is an 
<A HREF="http://www.boost.org/libs/graph/doc/using_property_maps.html">interior property</A>, 
that is, a
<A HREF="http://www.boost.org/libs/property_map/doc/property_map.html">property map</A> 
for a `vertex_point` can be extracted from any model of a `HalfedgeGraph`
using the <span class="textsc">Bgl</span>
<A HREF="http://www.boost.org/libs/graph/doc/PropertyGraph.html">PropertyGraph</A> interface:
`boost::get(vertex_point,graph)`

A point of a vertex can be directly accessed via:
- `boost::get(vertex_point,graph,vertex)`
- `boost::put(vertex_point,graph,vertex,newpoint)`
*/
enum vertex_point_t { vertex_point };

}

namespace boost {

/*!
\ingroup PkgBGLEnums
The constant `edge_index` is a 
<A HREF="http://www.boost.org/libs/graph/doc/PropertyTag.html">property tag</A> which identifies the <I>index</I> property
of an edge of a \sc{Bgl} 
<A HREF="http://www.boost.org/libs/graph/doc/Graph.html">Graph</A>.

*/
enum edge_index_t { edge_index };
}
