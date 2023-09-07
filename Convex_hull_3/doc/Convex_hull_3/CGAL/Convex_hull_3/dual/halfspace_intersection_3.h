namespace CGAL {

/*!
\ingroup PkgConvexHull3Functions

\brief computes robustly the intersection of the halfspaces defined by the planes contained in the range [`begin`, `end`) without constructing the dual points. The result is stored in the polyhedron `pm`.
If `origin` is given then it must be a point strictly inside the polyhedron. If an interior point is not given then it is computed using the function `halfspace_intersection_interior_point_3()` based on solving a linear program and thus is slower.

This version does not construct the dual points explicitly but uses a special traits class for the function `CGAL::convex_hull_3()` to handle predicates on dual points without constructing them.

Halfspaces are considered as lower halfspaces, that is if the plane equation is \f$ a\, x +b\, y +c\, z + d = 0 \f$ then the corresponding halfspace is defined by \f$ a\, x +b\, y +c\, z + d \le 0 \f$ .


\pre The point type of `origin` and the point type of the vertices of `PolygonMesh` must come from the same \cgal %Kernel.\pre if provided, `origin` is inside the intersection of halfspaces defined by the range `[begin, end)`.
\pre The computed intersection must be a bounded convex polyhedron.


\tparam PlaneIterator must be an input iterator where the value type is a model of the
        concept `Kernel::Plane_3` and this plane type must come from the same kernel as the point type.

\tparam PolygonMesh must be a model of `MutableFaceGraph`.

\sa `halfspace_intersection_with_constructions_3()`
 */

template <class PlaneIterator, class PolygonMesh>
void halfspace_intersection_3 (PlaneIterator begin, PlaneIterator end,
                               PolygonMesh &pm,
                               std::optional<Kernel_traits<std::iterator_traits<PlaneIterator>::value_type>::Kernel::Point_3> > origin = std::nullopt);

} /* namespace CGAL */
