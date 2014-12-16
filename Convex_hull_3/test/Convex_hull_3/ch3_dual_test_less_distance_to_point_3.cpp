#include <cassert>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Convex_hull_3/dual/Convex_hull_traits_dual_3.h>
#include <CGAL/Convex_hull_traits_3.h>

#include "include/to_dual.h"

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Plane_3 Plane;

typedef CGAL::Convex_hull_3::Convex_hull_traits_dual_3<K> Hull_traits_dual;
typedef Hull_traits_dual::Less_distance_to_point_3 Less_distance_to_point_3_dual;

typedef CGAL::Convex_hull_traits_3<K> Hull_traits;
typedef Hull_traits::Less_distance_to_point_3 Less_distance_to_point_3;

int main () {
    Hull_traits_dual traits_dual;
    Hull_traits traits;

    Less_distance_to_point_3_dual less_distance_to_point_dual = traits_dual.less_distance_to_point_3_object();
    Less_distance_to_point_3 less_distance_to_point = traits.less_distance_to_point_3_object();

    Plane p(1, 0, 0, -1), q(3, 0, 0, -1);
    Plane r(2, 0, 0, -1), rr(4, 0, 0, -1);

    // True
    assert(less_distance_to_point_dual(p, q, rr) == true);
    assert(less_distance_to_point_dual(p, q, rr) == less_distance_to_point(to_dual<K>(p), to_dual<K>(q), to_dual<K>(rr)));

    // False
    assert(less_distance_to_point_dual(p, q, r) == false);
    assert(less_distance_to_point_dual(p, q, r) == less_distance_to_point(to_dual<K>(p), to_dual<K>(q), to_dual<K>(r)));

    // Tests if r is one of the two first points
    // False
    assert(less_distance_to_point_dual(p, q, p) == false);
    assert(less_distance_to_point_dual(p, q, p) == less_distance_to_point(to_dual<K>(p), to_dual<K>(q), to_dual<K>(p)));

    assert(less_distance_to_point_dual(p, q, q) == false);
    assert(less_distance_to_point_dual(p, q, q) == less_distance_to_point(to_dual<K>(p), to_dual<K>(q), to_dual<K>(q)));

    return 0;
}

