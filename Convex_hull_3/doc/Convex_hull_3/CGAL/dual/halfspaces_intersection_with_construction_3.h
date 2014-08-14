namespace CGAL {

/*!
\ingroup PkgConvexHull3Functions

\brief computes the intersection of the halfspaces defined by the planes contained in the range [`begin`, `end`). The result is stored in the polyhedron `P`.
In order to do that, it is necessary to give the function a point inside the polyhedron named `origin` which is `CGAL::ORIGIN` by default.
This version constructs explicitly the dual points.

\attention Halfspaces are considered as lower halfspaces that is to say if the plane's equation is \f$ a\, x +b\, y +c\, z + d = 0 \f$ then the corresponding halfspace is defined by \f$ a\, x +b\, y +c\, z + d \le 0 \f$ .

\pre `origin` is inside the intersection of halfspaces defined by the range [`begin`, `end`).
\pre The intersection exists that is to say that it is a polyhedron.

\tparam PlaneIterator must be an input iterator with a value type  equivalent to `Polyhedron::Traits::Plane_3`.
\tparam Polyhedron must be a model of `ConvexHullPolyhedron_3`.
 */

template <class PlaneIterator, class Polyhedron>
void halfspaces_intersection_with_construction_3(PlaneIterator pbegin,
                                                 PlaneIterator pend,
                                                 Polyhedron &P,
                                                 typename Polyhedron::Traits::Point_3 const& origin = typename Polyhedron::Traits::Point_3(CGAL::ORIGIN));

} /* namespace CGAL */
