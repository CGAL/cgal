#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Voronoi_covariance_3/Convex_hull_traits_dual_3.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Plane_3 Plane;

typedef CGAL::Voronoi_covariance_3::Convex_hull_traits_dual_3<K, true> Hull_traits_dual;
typedef Hull_traits_dual::Equal_3 Equal_3;

int main (void) {
    Hull_traits_dual dual_traits;

    Equal_3 equal = dual_traits.equal_3_object();

    Plane p(1, 1, 1, 1), q(1, 1, 1, 1);
    std::cout << equal(p, q) << std::endl;
    assert(equal(p, q) == true);

    return 0;
}

