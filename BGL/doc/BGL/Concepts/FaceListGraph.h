/*!
\ingroup PkgBGLConcepts
\cgalConcept

The concept `FaceListGraph` refines the concept `FaceGraph` and adds
the requirement for traversal of all faces in a graph.

\cgalRefines `FaceGraph`
\cgalHasModel `CGAL::Polyhedron_3`
\cgalHasModel `CGAL::Surface_mesh`

\cgalHeading{Notations}

<dl>
<dt>`G`</dt> 	<dd>A type that is a model of `FaceListGraph`.</dd>
<dt>`g`</dt> 	<dd>An object of type `G`.</dd>
</dl>

\cgalHeading{Associated Types}

Type              | Description
----------------- | -----------
`boost::graph_traits<G>::%face_iterator`   | %Iterator over all faces.
`boost::graph_traits<G>::%faces_size_type` | Unsigned integer type for number of faces.


\cgalHeading{Valid Expressions}

Expression        |  returns                               | Description
----------------- | ---------------                        | -----------------------
`faces(g)`        |  `std::pair<face_iterator, face_iterator>` | An iterator range over all faces. 
`num_faces(g)`    |  `faces_size_type`                     | An upper bound of the number of faces of the graph.

\attention `num_faces()` may return a number larger than `std::distance(faces(g).first, faces(g).second)`.
This is the case for implementations only marking faces deleted in the face container.
<!--
This is for example the case for `CGAL::Surface_mesh` or `OpenMesh::PolyMesh_ArrayKernelT`. 
-->

*/
class FaceListGraph{};
