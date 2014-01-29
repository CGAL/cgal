
/*!
\ingroup PkgSurfaceMeshSimplificationConcepts
\cgalConcept

The concept `EdgeCollapsableMesh` describes the requirements for the type of 
triangulated surface mesh that can be passed to the 
simplification algorithm. 

The surface must be structurally equivalent to a polyhedral surface 
having only triangular faces. 
It can have any number of connected components, boundaries 
(borders and holes) and handles (arbitrary genus). 

\cgalRefines `HalfedgeGraph` 

\cgalHeading{Valid Expressions}

Let `v0v1` an edge of the triangulated surface mesh `ecm` and
`v0` and `v1` being the source and target vertices of that edge.
The mesh simplification algorithm requires the call to the function `halfedge_collapse(v0v1,ecm)`
to be valid and to return the vertex not removed after collapsing
the undirected edge `(v0v1,v1v0)`.

For `e` \f$ \in \{\f$ `v0v1,v1v0` \f$ \}\f$, let `en` and `ep` be the next and previous 
edges, that is `en = next_edge(e, mesh)`, `ep = prev_edge(e,mesh)`, and let 
`eno` and `epo` be their opposite edges, that is 
`eno = opposite_edge(en, mesh)` and `epo = opposite_edge(ep,mesh)`. 

Then, after the collapse of `(v0v1,v1v0)` the following holds: 

<UL> 
<LI>The edge `e` is no longer in `mesh`. 
<LI>One of \f$ \{\f$`v0,v1`\f$ \}\f$ is no longer in `mesh` while the other remains. 
\cgalFootnote{Most of the time v0 is the vertex being removed but in some cases removing the edge e requires v1 to be removed. See Figure \ref CollapseFigure5.}
Let `vgone` be the removed vertex and `vkept` be the remaining vertex. 
<LI>If `e` was a border edge, that is `get(is_border, e, mesh) == true`, then `next_edge(ep) == en`, and `prev_edge(en) == ep`. 
<LI>If `e` was not a border edge, that is `get(is_border, e, mesh) == false`, then `ep` and `epo` are no longer in `mesh` while `en` and `eno` are kept in `mesh`. 
<LI>For all edges `ie` in `in_edges(vgone,mesh)`, `target(ie,mesh) == vkept` and `source(opposite_edge(ie),mesh) == vkept`. 
<LI>No other incidence information has changed in `mesh`. 
</UL>


\image html general_collapse.png
\image latex general_collapse.png
<center><b>
General case. The following mesh elements are removed: triangles (\f$ v0,v1,vL\f$) and (\f$ v1,v0,vR\f$), edges \f$ (e,e')\f$, \f$ (ep,epo)\f$ and \f$ (ep',epo')\f$, and vertex \f$ v0\f$. 
</b></center>

\image html border_collapse3.png "When the collapsing edge is not itself a border, but is incident upon a border edge that is removed, the operation is the same as in the general case."
\image latex border_collapse3.png "When the collapsing edge is not itself a border, but is incident upon a border edge that is removed, the operation is the same as in the general case."

\image html border_collapse2.png
\image latex border_collapse2.png
<center><b>
When the collapsing edge is not itself a border, but is incident upon
a border edge that is <I>not</I> removed, the operation is still the
same as in the general case.
</b></center>

\image html border_collapse1.png
\image latex border_collapse1.png
<center><b>
When the collapsing edge is itself a border, only 1 triangle is
removed. Thus, even if \f$ (ep',epo')\f$ exists, it's not removed.
</b></center>

\anchor CollapseFigure5
\image html border_collapse4.png
\image latex border_collapse4.png
<center><b>
This figure illustrates the single exceptional case when removing \f$
(v0,v1)\f$ neccesarily implies removing \f$ (v1)\f$, thus \f$ (v0)\f$
remains.
</b></center>

\cgalHasModel `CGAL::Polyhedron_3<Traits>` (If it has only triangular faces),
using the specialization \link BGLPolyGT `boost::graph_traits< CGAL::Polyhedron_3<Traits> >` \endlink.

\sa \link BGLPolyGT `boost::graph_traits< CGAL::Polyhedron_3<Traits> >` \endlink
\sa `CGAL::halfedge_graph_traits< CGAL::Polyhedron_3<Traits> >`

*/

class EdgeCollapsableMesh {
public:
}; /* end EdgeCollapsableMesh */

/*!
Collapses the undirected edge `(v0v1,v1v0)` replacing it with `v0` or `v1`, 
as described in the paragraph above.
\pre This function requires `mesh` to be an oriented 2-manifold with or without boundaries. Furthermore, the undirected edge `(v0v1,v1v0)` must satisfy the <I>link condition</I> \cgalCite{degn-tpec-98}, which guarantees that the surface is also 2-manifold after the edge collapse. 
\relates EdgeCollapsableMesh
*/ 
template<class EdgeCollapsableMesh> 
typename boost::graph_traits<EdgeCollapsableMesh>::vertex_descriptor 
halfedge_collapse(typename boost::graph_traits<EdgeCollapsableMesh>::edge_descriptor const& ue, EdgeCollapsableMesh& mesh); 


