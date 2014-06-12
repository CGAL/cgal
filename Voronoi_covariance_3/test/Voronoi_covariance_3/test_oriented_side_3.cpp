#include <cassert>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Voronoi_covariance_3/predicates.h>
#include <CGAL/Convex_hull_traits_3.h>

#include "include/to_dual.h"

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Plane_3 Plane;

typedef CGAL::Voronoi_covariance_3::Oriented_side_3_dual_point<K> Oriented_side_3_dual;
typedef CGAL::Voronoi_covariance_3::Plane_dual<K> Plane_dual;

typedef CGAL::Convex_hull_traits_3<K> Hull_traits;
typedef Hull_traits::Oriented_side_3 Oriented_side_3;

void printOrientation (CGAL::Oriented_side o) {
    if (o == CGAL::ON_POSITIVE_SIDE) {
        std::cout << "positive" << std::endl;
    } else if (o == CGAL::ON_NEGATIVE_SIDE) {
        std::cout << "negative" << std::endl;
    } else {
        std::cout << "boundary" << std::endl;
    }
}

int main (void) {
    Oriented_side_3_dual oriented_dual;
    Oriented_side_3 oriented;

    Plane a(1, 1, 0, -1), b(1, 2, 1, -1), c(1, 3, 0, -1); // dual plane : x = 1
    Plane_dual p(a,b,c);
    Plane q(0, 0, 0, -1), qq(2, 0, 0, -1), qqq(1, 1, 0, -1);

    // p1.d() * q.d() > 0
    // Positive
    /* printOrientation(oriented_dual(p, q)); */
    /* printOrientation(oriented(to_dual_plane<K>(p), to_dual<K>(q))); */
    assert(oriented_dual(p, qq) == CGAL::ON_POSITIVE_SIDE);
    assert(oriented_dual(p, q) == oriented(to_dual_plane<K>(p), to_dual<K>(q)));

    // Negative
    /* printOrientation(oriented_dual(p, qq)); */
    /* printOrientation(oriented(to_dual_plane<K>(p), to_dual<K>(qq))); */
    assert(oriented_dual(p, q) == CGAL::ON_NEGATIVE_SIDE);
    assert(oriented_dual(p, qq) == oriented(to_dual_plane<K>(p), to_dual<K>(qq)));

    // Boundary
    /* printOrientation(oriented_dual(p, qqq)); */
    /* printOrientation(oriented(to_dual_plane<K>(p), to_dual<K>(qqq))); */
    assert(oriented_dual(p, qqq) == CGAL::ON_oriented_dual_BOUNDARY);
    assert(oriented_dual(p, qqq) == oriented(to_dual_plane<K>(p), to_dual<K>(qqq)));

    // p1.d() * q.d() < 0
    Plane aa(-1, -1, 0, 1);
    Plane_dual pp(aa, b, c);

    // Positive
    /* printOrientation(oriented_dual(pp, q)); */
    /* printOrientation(oriented(to_dual_plane<K>(pp), to_dual<K>(q))); */
    assert(oriented_dual(pp, qq) == CGAL::ON_POSITIVE_SIDE);
    assert(oriented_dual(pp, q) == oriented(to_dual_plane<K>(pp), to_dual<K>(q)));

    // Negative
    /* printOrientation(oriented_dual(pp, qq)); */
    /* printOrientation(oriented(to_dual_plane<K>(pp), to_dual<K>(qq))); */
    assert(oriented_dual(pp, q) == CGAL::ON_NEGATIVE_SIDE);
    assert(oriented_dual(pp, qq) == oriented(to_dual_plane<K>(pp), to_dual<K>(qq)));

    // Boundary
    /* printOrientation(oriented_dual(pp, qqq)); */
    /* printOrientation(oriented(to_dual_plane<K>(pp), to_dual<K>(qqq))); */
    assert(oriented_dual(pp, qqq) == CGAL::ON_oriented_dual_BOUNDARY);
    assert(oriented_dual(pp, qqq) == oriented(to_dual_plane<K>(pp), to_dual<K>(qqq)));

    return 0;
}

