namespace CGAL {

/*!
\ingroup PkgTriangulation3

fills a face graph with the link of a triangulation vertex. 


\pre `T.dimension()`==3.

\tparam Triangulation must be a \cgal 3D triangulation. 
\tparam TriangleMesh must be a model of the concept `MutableFaceGraph`. 

\param t the 3D triangulation
\param vh the vertex handle of the vertex of the link
\param tm the triangle mesh
\param no_infinite_faces If `vh` is on the convex hull
         of the triangulation, `no_infinite_faces == true` generates a triangle mesh with a border.

\returns the vertex descriptor of the `tm` corresponding to the infinite vertex of `t`, 
         if `vh` is on the convex hull of the triangulation, and if `no_infinite_faces == false`.

\sa `convex_hull_3_to_polyhedron_3()`
 
*/
template <class Triangulation, class TriangleMesh>
typename boost::graph_trait<FG>::vertex_descriptor
link_to_face_graph(const Triangulation& t,
                   typename Triangulation_3::Vertex_handle vh,
                   TriangleMesh& tm, 
                   bool no_infinite_faces = true);

} /* namespace CGAL */
