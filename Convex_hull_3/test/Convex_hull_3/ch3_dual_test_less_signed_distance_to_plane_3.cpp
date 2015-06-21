#include <cassert>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Convex_hull_3/dual/Convex_hull_traits_dual_3.h>
#include <CGAL/Convex_hull_traits_3.h>

#include "include/to_dual.h"

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Plane_3 Plane;
typedef CGAL::Convex_hull_3::Plane_dual<K> Plane_dual;

typedef CGAL::Convex_hull_3::Convex_hull_traits_dual_3<K> Hull_traits_dual;
typedef Hull_traits_dual::Less_signed_distance_to_plane_3 Less_signed_distance_to_plane_3_dual;

typedef CGAL::Convex_hull_traits_3<K> Hull_traits;
typedef Hull_traits::Less_signed_distance_to_plane_3 Less_signed_distance_to_plane_3;

int main () {
    Hull_traits_dual traits_dual;
    Hull_traits traits;

    Less_signed_distance_to_plane_3 less_signed_distance = traits.less_signed_distance_to_plane_3_object();
    Less_signed_distance_to_plane_3_dual less_signed_distance_dual = traits_dual.less_signed_distance_to_plane_3_object();

    // dual plane : x = 1
    Plane a(1, 1, 0, -1), b(1, 2, 1, -1), c(1, 3, 0, -1);
    Plane_dual p(a,b,c);

    Plane q(2, 0, 0, -1), r(3, 0, 0, -1) , rr(0.5, 0, 0, -1);

    // False
    assert(less_signed_distance_dual(p, q, r) == false);
    assert(less_signed_distance_dual(p, q, r) == less_signed_distance(to_dual_plane<K>(p), to_dual<K>(q), to_dual<K>(r)));

    // True
    assert(less_signed_distance_dual(p, q, rr) == true);
    assert(less_signed_distance_dual(p, q, rr) == less_signed_distance(to_dual_plane<K>(p), to_dual<K>(q), to_dual<K>(rr)));

    return 0;
}

