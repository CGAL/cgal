#ifndef TO_DUAL_H
#define TO_DUAL_H

// Convert a dual point to a point
template <typename R>
typename R::Point_3 to_dual (typename R::Plane_3 const& p) {
    typename R::Point_3 pp(-p.a() / p.d(), -p.b() / p.d(), -p.c() / p.d());

    return pp;
}

// Convert a dual plane to a plane
template <typename R>
CGAL::Point_triple<R> to_dual_plane (CGAL::Convex_hull_3::Plane_dual<R> const& p) {
    typename R::Plane_3 p1 = p.p1;
    typename R::Plane_3 p2 = p.p2;
    typename R::Plane_3 p3 = p.p3;

    CGAL::Point_triple<R> pp(to_dual<R>(p1),
                             to_dual<R>(p2),
                             to_dual<R>(p3));

    return pp;
}

#endif

