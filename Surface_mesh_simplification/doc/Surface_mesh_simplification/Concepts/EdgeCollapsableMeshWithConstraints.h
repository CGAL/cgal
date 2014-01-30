
/*!
\ingroup PkgSurfaceMeshSimplificationConcepts
\cgalConcept

The concept `EdgeCollapsableMeshWithConstraints` describes additional requirements
for the type of triangulated surface meshes that can be passed to the
simplification algorithm.

\cgalRefines `EdgeCollapsableMesh`

\cgalHeading{Valid Expressions}

Let `v0v1` be an edge of the triangulated surface mesh `ecm` and
`v0` and `v1` being the source and target vertices of that edge.
The mesh simplification algorithm requires the call to the function `halfedge_collapse(v0v1,ecm)`
to be valid and to return the vertex not removed after collapsing
the undirected edge `(v0v1,v1v0)`.

For `e` \f$ \in \{\f$ `v0v1,v1v0` \f$ \}\f$, let `en` and `ep` be the next and previous
edges, that is `en = next_edge(e, mesh)`, `ep = prev_edge(e,mesh)`, and let
`eno` and `epo` be their opposite edges, that is
`eno = opposite_edge(en, mesh)` and `epo = opposite_edge(ep,mesh)`.

Then, after the collapse of `(v0v1,v1v0)` the invariants described in the concept `EdgeCollapsableMesh` hold
if `ep` is not constrained. Otherwise, it is `en` that is removed from `ecm`.

\image html collapse_constraints.png
\image latex collapse_constraints.png

\cgalHasModel `CGAL::Polyhedron_3<Traits>` (If it has only triangular faces),
using the specialization \link BGLPolyGT `boost::graph_traits< CGAL::Polyhedron_3<Traits> >` \endlink.

\sa \link BGLPolyGT `boost::graph_traits< CGAL::Polyhedron_3<Traits> >` \endlink
\sa `CGAL::halfedge_graph_traits< CGAL::Polyhedron_3<Traits> >`

*/

class EdgeCollapsableMeshWithConstraints {
public:
}; /* end EdgeCollapsableMeshWithConstraints */

/*!
Collapses the undirected edge `(v0v1,v1v0)` replacing it with `v0` or `v1`,
as described in the paragraph above and guarantees that an halfedge `he`, for which `get(edge_is_constrained_map, he)==true`, is not removed after the collapse.
\tparam EdgeCollapsableMeshWithConstraints a model of `HalfedgeGraph`
\tparam EdgeIsConstrainedMap a model of `ReadablePropertyMap` with the edge descriptor of
       `EdgeCollapsableMeshWithConstraints` as key type and a boolean as value type.
        It indicates if an edge is constrained or not.
\pre This function requires `mesh` to be an oriented 2-manifold with or without boundaries. Furthermore, the undirected edge `(v0v1,v1v0)` must satisfy the <I>link condition</I> \cgalCite{degn-tpec-98}, which guarantees that the surface is also 2-manifold after the edge collapse.
\pre `get(edge_is_constrained_map, v0v1)==get(edge_is_constrained_map, v1v0)==false`.
\pre `v0` and `v1` are not both incident to a constrained edge.
\relates EdgeCollapsableMeshWithConstraints
*/
template<class EdgeCollapsableMeshWithConstraints,class EdgeIsConstrainedMap>
typename boost::graph_traits<EdgeCollapsableMeshWithConstraints>::vertex_descriptor
halfedge_collapse(typename boost::graph_traits<EdgeCollapsableMeshWithConstraints>::edge_descriptor const& ue,
                  EdgeCollapsableMeshWithConstraints& mesh,
                  EdgeIsConstrainedMap edge_is_constrained_map);


