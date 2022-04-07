#include <cassert>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Convex_hull_3/dual/Convex_hull_traits_dual_3.h>
#include <CGAL/Convex_hull_traits_3.h>

#include "include/to_dual.h"

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Plane_3 Plane;

typedef CGAL::Convex_hull_3::Convex_hull_traits_dual_3<K> Hull_traits_dual;
typedef Hull_traits_dual::Coplanar_3 Coplanar_3_dual;

typedef CGAL::Convex_hull_traits_3<K> Hull_traits;
typedef Hull_traits::Coplanar_3 Coplanar_3;

int main () {
    Hull_traits_dual traits_dual;
    Hull_traits traits;

    Coplanar_3_dual coplanar_dual = traits_dual.coplanar_3_object();
    Coplanar_3 coplanar = traits.coplanar_3_object();

    // Simple true test with equal planes
    K::Plane_3 p(1, 1, 1, 1), q(1, 1, 1, 1), r(1, 1, 1, 1), s(1, 1, 1, 1);
    assert(coplanar_dual(p, q, r, s) == true);
    assert(coplanar_dual(p, q, r, s) == coplanar(to_dual<K>(p), to_dual<K>(q), to_dual<K>(r), to_dual<K>(s)));

    // Simple true test
    K::Plane_3 p2(1, 1, 1, -1), q2(2, 1, 1, -1), r2(2, 2, 1, -1), s2(1, 2, 1, -1);
    assert(coplanar_dual(p2, q2, r2, s2) == true);
    assert(coplanar_dual(p2, q2, r2, s2) == coplanar(to_dual<K>(p2), to_dual<K>(q2), to_dual<K>(r2), to_dual<K>(s2)));

    // Simple false test
    K::Plane_3 p3(1, 1, 1, -1), q3(2, 1, 1, -1), r3(2, 2, 1, -1), s3(1, 1, 2, -1);
    assert(coplanar_dual(p3, q3, r3, s3) == false);
    assert(coplanar_dual(p3, q3, r3, s3) == coplanar(to_dual<K>(p3), to_dual<K>(q3), to_dual<K>(r3), to_dual<K>(s3)));

    return 0;
}

