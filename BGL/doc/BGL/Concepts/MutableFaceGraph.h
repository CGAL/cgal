/*!
\ingroup PkgBGLConcepts
\cgalConcept

The concept `MutableFaceGraph` refines the concepts `FaceGraph` and `MutableHalfedgeGraph` and adds
the requirement for operations to add faces and to modify face-halfedge relations.

\cgalRefines `FaceGraph`
\cgalRefines `MutableHalfedgeGraph`
\cgalHasModel `CGAL::Polyhedron_3`
\cgalHasModel `CGAL::Surface_mesh`

\cgalHeading{Notations}

<dl>
<dt>`G`</dt> 	<dd>A type that is a model of `MutableFaceGraph`.</dd>
<dt>`g`</dt> 	<dd>An object of type `G`.</dd>
<dt>`h`</dt> 	<dd>A halfedge descriptor.</dd>
<dt>`f`</dt> 	<dd>A face descriptor.</dd>
</dl>

\cgalHeading{Valid Expressions}

Expression              | returns           | Description                           
----------------------- | ------------      | -----------
`add_face(g)`           | `face_descriptor` | Adds a new face to the graph with no corresponding halfedge set.
`remove_face(f, g)`     | `void`            | Removes `f` from the graph.
`set_face(h, f, g)`     | `void`            | Sets the corresponding face of `h` to `f`.
`set_halfedge(f, h, g)` | `void`            | Sets the corresponding halfedge of `f` to `h`.

*/
class MutableFaceGraph{};
