/*!
\ingroup PkgBGLConcepts
\cgalConcept

The concept `HalfedgeListGraph` refines the concept `HalfedgeGraph`
and adds the requirements for traversal of all halfedges in the graph.

\cgalRefines `HalfedgeGraph`
\cgalHasModel `CGAL::Polyhedron_3`
\cgalHasModel `CGAL::Surface_mesh`

\cgalHeading{Notations}

<dl>
<dt>`G`</dt> 	<dd>A type that is a model of `HalfedgeListGraph`.</dd>
<dt>`g`</dt> 	<dd>An object of type `G`.</dd>
</dl>

\cgalHeading{Associated Types}

Type                 | Description
-------------------- | ------------
`boost::graph_traits<G>::%halfedge_iterator`      | A `BidirectionalIterator` over all halfedges in a graph. Must be `DefaultConstructible`, `Assignable`, `EqualityComparable`.
`boost::graph_traits<G>::%halfedges_size_type`    | A size type.


\cgalHeading{Valid Expressions}

Expression                            | Returns                                   | Description  
------------------------------------- | ------------------------------------------| -----------
`num_halfedges(g)`                    | `halfedges_size_type`                     | An upper bound of the number of halfedges of the graph.
`halfedges(g)`                        | `std::pair<halfedge_iterator,halfedge_iterator>` | An iterator range over the halfedges of the graph.

\attention `num_halfedges()` may return a number larger than `std::distance(halfedges(g).first,halfedges(g).second)`.
This is the case for implementations only marking halfedges deleted in the halfedge container.

<!--
This is for example the case for `CGAL::Surface_mesh` or `OpenMesh::PolyMesh_ArrayKernelT`. 
-->

*/

class HalfedgeListGraph {};
