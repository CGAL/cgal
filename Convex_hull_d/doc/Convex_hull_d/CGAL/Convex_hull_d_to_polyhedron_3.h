namespace CGAL {

/*!
\deprecated This package is deprecated since the version 4.6 of \cgal. The package \ref PkgTriangulations should be used instead.

converts the convex hull `C` to polyhedral surface stored in
`P`.
\pre `dim == 3` and `dcur == 3`.

\relates CGAL::Convex_hull_d
*/
template <class R, class T, class HDS>
void convex_hull_d_to_polyhedron_3( const Convex_hull_d<R>& C, Polyhedron_3<T,HDS>& P) ;

/*!
\deprecated This package is deprecated since the version 4.6 of \cgal. The package \ref PkgTriangulations should be used instead

constructs the representation of the surface of `C` as a
bidirected LEDA graph `G`.
\pre `dim == 3`.

\relates CGAL::Convex_hull_d
*/
template <class R>
void d3_surface_map(const Convex_hull_d<R>& C, GRAPH< typename Convex_hull_d<R>::Point_d ,int>& G);

} /* end namespace CGAL */
