/*!
\ingroup PkgBGLConcepts
\cgalConcept

The concept `FaceGraph` refines the concept `HalfedgeGraph`. 
It adds the requirements for a graph to explicitly
maintain faces described by halfedges, to provide access from a face to
an incident halfedge, and to provide access from a halfedge to its incident
face. 

\cgalRefines `HalfedgeGraph`
\cgalHasModel `CGAL::Polyhedron_3`
\cgalHasModel `CGAL::Surface_mesh`

\cgalHeading{Notations}

<dl>
<dt>`G`</dt> 	<dd>A type that is a model of `FaceGraph`.</dd>
<dt>`g`</dt> 	<dd>An object of type `G`.</dd>
<dt>`e`</dt> 	<dd>An edge descriptor.</dd>
<dt>`f`</dt> 	<dd>A face descriptor.</dd>
<dt>`h`</dt> 	<dd>A halfedge descriptor.</dd>
</dl>

\cgalHeading{Associated Types}

Type                 | Description
-------------------- | ------------
`boost::graph_traits<G>::%face_descriptor`              | A `face_descriptor` corresponds to a unique face in a graph. Must be `DefaultConstructible`, `Assignable`, `EqualityComparable` and `LessThanComparable`.


\cgalHeading{Valid Expressions}

Expression                             | Returns                                                                  | Description  
-------------------------------------- | ------------------------------------------------------------------------ | ------------------------
`face(h, g)`                           | `face_descriptor`                                                        | The face incident to halfedge `h`.
`halfedge(f, g)`                       | `halfedge_descriptor`                                                    | A halfedge incident to face `f`.
`degree(f,g)`                          | `degree_size_type`                                                       | The number of halfedges, incident to face `f`.
`boost::graph_traits<G>::%null_face()` | `face_descriptor`                                                        | Returns a special face that is not equal to any other face.


*/
class FaceGraph {};
