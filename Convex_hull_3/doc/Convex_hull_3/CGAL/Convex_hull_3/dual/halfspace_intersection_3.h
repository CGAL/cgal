namespace CGAL {

/*!
\ingroup PkgConvexHull3Functions

\brief computes robustly the intersection of the halfspaces defined by the planes contained in the range [`begin`, `end`) without constructing the dual points. The result is stored in the polyhedron `P`.
If `origin` is given then it must be a point strictly inside the polyhedron. If an interior point is not given then it is computed using a linear program and thus is slower.

This version does not construct the dual points explicitely but uses a special traits class for the function `CGAL::convex_hull_3()` to handle predicates on dual points without constructing them.

\attention Halfspaces are considered as lower halfspaces that is to say if the plane's equation is \f$ a\, x +b\, y +c\, z + d = 0 \f$ then the corresponding halfspace is defined by \f$ a\, x +b\, y +c\, z + d \le 0 \f$ .
\attention The value type of `PlaneIterator` and the point type of `origin` must come from the same \cgal %Kernel.

\pre if provided, `origin` is inside the intersection of halfspaces defined by the range `[begin, end)`.
\pre The computed intersection must be a bounded convex polyhedron.

\tparam PlaneIterator must be an input iterator where the value type must be `Polyhedron::Traits::Plane_3`
\tparam Polyhedron must be a model of `ConvexHullPolyhedron_3`.

\sa `halfspace_intersection_with_constructions_3()` 
 */

template <class PlaneIterator, class Polyhedron>
void halfspace_intersection_3 (PlaneIterator begin, PlaneIterator end,
                               Polyhedron &P,
                               boost::optional<Polyhedron::Vertex::Point_3> origin = boost::none);

} /* namespace CGAL */
