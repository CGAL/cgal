/*!
\ingroup PkgBGLConcepts
\cgalConcept

The concept `HalfedgeGraph` is a refinement of the \sc{Bgl} concept
<A HREF="http://www.boost.org/libs/graph/doc/Graph.html">`Graph`</A> and adds the notion of a *halfedge*: Each edge is
associated with two *opposite* halfedges with source and target vertices swapped.
Furthermore, halfedges have a *successor* and *predecessor*,
and form cycles we call *faces*. However, this concept 
does not introduce a face type. 
A `HalfedgeGraph` is undirected and does not allow parallel edges.

Using the composition of the *successor* and *opposite* functions results 
in another cycle, namely the cycle of halfedges which are incident to
the same vertex. We refer to \ref PkgBGLIterators for a description of
iterators and circulators for these halfedge cycles.


\cgalRefines <A HREF="http://www.boost.org/libs/graph/doc/Graph.html">`Graph`</A>
\cgalRefines <A HREF="http://www.boost.org/libs/graph/doc/PropertyGraph.html">`PropertyGraph`</A>

A model of `HalfedgeGraph` must have the interior property `vertex_point` attached to its vertices.

\cgalHasModel `CGAL::Polyhedron_3`
\cgalHasModel `CGAL::Surface_mesh`

\cgalHeading{Notations}

<dl>
<dt>`G`</dt> 	  <dd>A type that is a model of `HalfedgeGraph`.</dd>
<dt>`g`</dt> 	  <dd>An object of type `G`.</dd>
<dt>`u`, `v`</dt> <dd>Vertex descriptors.</dd>
<dt>`e`</dt> 	  <dd>An edge descriptor.</dd>
<dt>`h`</dt> 	  <dd>A halfedge descriptor.</dd>
</dl>

\cgalHeading{Associated Types}

Type                                                      | Description
--------------------------------------------------------- | ------------
`boost::graph_traits<G>::%halfedge_descriptor`             | A `halfedge_descriptor` corresponds to a halfedge in a graph. Must be `DefaultConstructible`, `Assignable`, `EqualityComparable` and `LessThanComparable`.


\cgalHeading{Valid Expressions}

Expression                              | Returns                                                                      | Description  
--------------------------------------- | ---------------------------------------------------------------------------- | -----------
`edge(h, g)`                            | `edge_descriptor`                                                            | The edge corresponding to `h` and `opposite(h)`.
`halfedge(e, g)`                        | `halfedge_descriptor`                                                        | One of the halfedges corresponding to `e`.
`halfedge(v, g)`                        | `halfedge_descriptor`                                                        | A halfedge with target `v`. 
`halfedge(u, v, g)`                     | `std::pair<halfedge_descriptor,bool>`                                        | The halfedge with source `u` and target `v`. The Boolean is `true`, iff this halfedge exists.
`opposite(h, g)`                        | `halfedge_descriptor`                                                        | The halfedge with source and target swapped.
`source(h,g)`                           | `vertex_descriptor`                                                          | The source vertex of `h`.
`target(h,g)`                           | `vertex_descriptor`                                                          | The target vertex of `h`.
`next(h, g)`                            | `halfedge_descriptor`                                                        | The next halfedge around its face.
`prev(h, g)`                            | `halfedge_descriptor`                                                        | The previous halfedge around its face.
*/
class HalfedgeGraph {};
