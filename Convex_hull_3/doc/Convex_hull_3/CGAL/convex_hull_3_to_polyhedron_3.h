namespace CGAL {

/*!
\ingroup PkgConvexHull3Functions

fills a polyhedron with the convex hull of a set of 3D points contained in a 3D triangulation of \cgal. 

The polyhedron `P` is cleared and the convex hull of the set of 3D points is stored in `P`.

\attention This function does not compute the plane equations of the faces of `P`.

\pre `T.dimension()`==3.

\tparam Triangulation is a \cgal 3D triangulation. 
\tparam Polyhedron is an instantiation of `CGAL::Polyhedron_3<Traits>`. 

\sa `convex_hull_3` 
*/
template <class Triangulation, class Polyhedron>
void convex_hull_3_to_polyhedron_3(const Triangulation& T,Polyhedron& P);

} /* namespace CGAL */
