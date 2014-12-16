#include <cassert>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Convex_hull_3/dual/Convex_hull_traits_dual_3.h>
#include <CGAL/Convex_hull_traits_3.h>

#include "include/to_dual.h"

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Plane_3 Plane;

typedef CGAL::Convex_hull_3::Convex_hull_traits_dual_3<K> Hull_traits_dual;
typedef Hull_traits_dual::Collinear_3 Collinear_3_dual;

typedef CGAL::Convex_hull_traits_3<K> Hull_traits;
typedef Hull_traits::Collinear_3 Collinear_3;

int main () {
    Hull_traits_dual traits_dual;
    Hull_traits traits;

    Collinear_3_dual collinear_dual = traits_dual.collinear_3_object();
    Collinear_3 collinear = traits.collinear_3_object();

    // Simple true test with equal planes
    Plane p(1, 1, 1, 1), q(1, 1, 1, 1), r(1, 1, 1, 1);
    assert(collinear_dual(p, q, r) == true);
    assert(collinear_dual(p, q, r) == collinear(to_dual<K>(p), to_dual<K>(q), to_dual<K>(r)));

    // Simple true test with non equal planes
    Plane p2(1, 1, 1, -1), q2(2, 2, 2, -1), r2(3, 3, 3, -1);
    assert(collinear_dual(p2, q2, r2) == true);
    assert(collinear_dual(p2, q2, r2) == collinear(to_dual<K>(p2), to_dual<K>(q2), to_dual<K>(r2)));

    // Simple false test
    Plane p3(1, 1, 1, -1), q3(2, 2, 2, -1), r3(3, 3, 2, -1);
    assert(collinear_dual(p3, q3, r3) == false);
    assert(collinear_dual(p3, q3, r3) == collinear(to_dual<K>(p3), to_dual<K>(q3), to_dual<K>(r3)));

    // Almost but not collinear
    double almost1 = 0.9999999999;
    Plane p4(0, 0, 0, -1), q4(1, 1, 1, -1);
    Plane r4(almost1, almost1, almost1+1e-15, -1);
    assert(collinear_dual(p4, q4, r4) == false);
    assert(collinear_dual(p4, q4, r4) == collinear(to_dual<K>(p4), to_dual<K>(q4), to_dual<K>(r4)));

    return 0;
}

