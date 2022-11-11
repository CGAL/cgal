namespace CGAL {

/*!
\ingroup PkgConvexHull3Functions

fills a polyhedron with the convex hull of a set of 3D points contained in a 3D triangulation of \cgal.

The polyhedron `pm` is cleared and the convex hull of the set of 3D points is stored in `pm`.

\pre `T.dimension()`==3.

\tparam Triangulation must be a \cgal 3D triangulation
\tparam PolygonMesh must be a model of the concept `MutableFaceGraph`

\sa `convex_hull_3()`
\sa `link_to_face_graph()`

*/
template <class Triangulation, class PolygonMesh>
void convex_hull_3_to_face_graph(const Triangulation& T,PolygonMesh& pm);

} /* namespace CGAL */
