namespace CGAL {

/*!
\ingroup PkgConvexHull3Functions

The function `convex_hull_3_to_polyhedron_3` fills a polyhedron with the convex hull 
of a set of 3D points contained in a 3D triangulation of \cgal. 

Requirements 
-------------- 

This function requires the following: 
<OL> 
<LI>`Triangulation_3` is a \cgal\ 3D triangulation. 
<LI>`Polyhedron_3` is an instantiation of `CGAL::Polyhedron_3<Traits>`. 
</OL> 

\sa `CGAL::convex_hull_3` 

The polyhedron `P` is cleared and the convex hull of the set of 3D points is stored in `P`.
The plane equations of each face are not computed.
\pre `T.dimension()`==3.

*/

template <class Triangulation_3, class Polyhedron_3>
void convex_hull_3_to_polyhedron_3(const Triangulation_3& T,Polyhedron_3& P);

} /* namespace CGAL */
