namespace CGAL {

/*!
\ingroup PkgTriangulation3Ref

fills the face graph `tm` with the <A HREF="https://en.wikipedia.org/wiki/Simplicial_complex#Closure.2C_star.2C_and_link">link</a> of triangulation vertex `vh`.


\pre `T.dimension()`==3.

\tparam Triangulation must be a \cgal 3D triangulation.
\tparam TriangleMesh must be a model of the concept `MutableFaceGraph`.

\param t the 3D triangulation
\param vh the vertex handle of the vertex
\param tm the triangle mesh
\param no_infinite_faces If `vh` is on the convex hull
         of the triangulation, `no_infinite_faces == true` generates a triangle mesh with a border.
         Otherwise, this parameter is ignored.

\returns the vertex descriptor of the triangle mesh `tm` corresponding to the infinite vertex of `t`,
         if `vh` is on the convex hull of the triangulation, and if `no_infinite_faces == false`.
         Otherwise, an arbitrary vertex descriptor of the triangle mesh `tm`.

\sa `convex_hull_3_to_face_graph()`

*/
template <class Triangulation, class TriangleMesh>
typename boost::graph_trait<FG>::vertex_descriptor
link_to_face_graph(const Triangulation& t,
                   typename Triangulation::Vertex_handle vh,
                   TriangleMesh& tm,
                   bool no_infinite_faces = true);

} /* namespace CGAL */
