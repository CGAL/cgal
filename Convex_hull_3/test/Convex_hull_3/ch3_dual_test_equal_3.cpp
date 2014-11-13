#include <cassert>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Convex_hull_3/dual/Convex_hull_traits_dual_3.h>
#include <CGAL/Convex_hull_traits_3.h>

#include "include/to_dual.h"

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Plane_3 Plane;
typedef K::Point_3 Point;

typedef CGAL::Convex_hull_3::Convex_hull_traits_dual_3<K> Hull_traits_dual;
typedef Hull_traits_dual::Equal_3 Equal_3_dual;

typedef CGAL::Convex_hull_traits_3<K> Hull_traits;
typedef Hull_traits::Equal_3 Equal_3;

int main (void) {
    Point origin(0, 0, 0);

    Hull_traits_dual traits_dual(origin);
    Hull_traits traits;

    Equal_3_dual equal_dual = traits_dual.equal_3_object();
    Equal_3 equal = traits.equal_3_object();

    // True test with equal planes
    Plane p(1, 1, 1, 1), q(1, 1, 1, 1);
    assert(equal_dual(p, q) == true);
    assert(equal_dual(p, q) == equal(to_dual<K>(p), to_dual<K>(q)));

    // True test with planes whose coefficients are multiple of the other
    Plane pp(1, 1, 1, 1), qq(2, 2, 2, 2);
    assert(equal_dual(pp, qq) == true);
    assert(equal_dual(pp, qq) == equal(to_dual<K>(pp), to_dual<K>(qq)));

    // False test with non equal planes
    Plane p2(1, 2, 1, 1), q2(2, 2, 2, 2);
    assert(equal_dual(p2, q2) == false);
    assert(equal_dual(p2, q2) == equal(to_dual<K>(p2), to_dual<K>(q2)));

    // True test using exact computation
    double x = 1.123456789, y = 1.987654321, e=1e-15;
    Plane ppp(x, x, x, x), qqq(y+e, y+e, y+e, y+e);
    assert(equal_dual(ppp, qqq) == true);
    assert(equal_dual(ppp, qqq) == equal(to_dual<K>(ppp), to_dual<K>(qqq)));

    return 0;
}

