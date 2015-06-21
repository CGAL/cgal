/*!
\ingroup PkgBGLConcepts
\cgalConcept

The concept `MutableHalfedgeGraph` refines the concept `HalfedgeGraph`
and adds the requirements for operations to add vertices and edges, and to
update the incidence information between vertices and halfedges.


\cgalRefines `HalfedgeGraph`
\cgalHasModel `CGAL::Polyhedron_3`
\cgalHasModel `CGAL::Surface_mesh`

\cgalHeading{Notations}

<dl>
<dt>`G`</dt>        <dd>A type that is a model of `MutableHalfedgeGraph`.</dd>
<dt>`g`</dt>        <dd>An object of type `G`.</dd>
<dt>`v`</dt>        <dd>A vertex descriptor.</dd>
<dt>`h`, `h1`, `h2`</dt> <dd>Halfedge descriptors.</dd>
<dt>`e`</dt>        <dd>An edge descriptor.</dd>
</dl>

\cgalHeading{Valid Expressions}

Expression                | returns             | Description                           
------------------------- | ------------        | -----------
`add_vertex(g)`           | `vertex_descriptor` | Adds a new vertex to the graph without initializing the connectivity.
`remove_vertex(v, g)`     | `void`              | Removes `v` from the graph.
`add_edge(g)`             | `edge_descriptor`   | Adds two opposite halfedges to the graph without initializing the connectivity.
`remove_edge(e, g)`       | `void`              | Removes the two halfedges corresponding to `e` from the graph.
`set_target(h, v, g)`     | `void`              | Sets the target vertex of `h` and the source of `opposite(h)` to `v`.
`set_halfedge(v, h, g)`   | `void`              | Sets the halfedge of `v` to `h`. The target vertex of `h` must be `v`. 
`set_next(h1, h2, g)`     | `void`              | Sets the successor of `h1` around a face to `h2`, and the prededecessor of `h2` to `h1`.


*/
class MutableHalfedgeGraph{};
