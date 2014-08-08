namespace CGAL {

/*!
\ingroup PkgConvexHull3Functions

\brief computes the intersection of the halfspaces defined by the planes contained in the range [`begin`, `end`). The result is stored in the polyhedron `P`.

\attention Halfspaces are considered as lower halfspaces that is to say if the plane's equation is \f$ a\, x +b\, y +c\, z + d = 0 \f$ then the corresponding halfspace is defined by \f$ a\, x +b\, y +c\, z + d \le 0 \f$ .

\pre `CGAL::ORIGIN` is inside the intersection of halfspaces defined by the range [`begin`, `end`)

\tparam PlaneIterator must be an input iterator with a value type  equivalent to `Polyhedron::Traits::Plane_3`.
\tparam Polyhedron must be a model of `ConvexHullPolyhedron_3`.
 */

template <class PlaneIterator, class Polyhedron>
void halfspaces_intersection_with_construction_3(PlaneIterator pbegin,
                                                 PlaneIterator pend,
                                                 Polyhedron &P);

} /* namespace CGAL */
