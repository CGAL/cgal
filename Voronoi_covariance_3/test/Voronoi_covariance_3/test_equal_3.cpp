#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Voronoi_covariance_3/predicates.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Plane_3 Plane;

typedef CGAL::Voronoi_covariance_3::Equal_3_dual_point<K> Equal_3_dual;

int main (void) {
    Equal_3_dual equal;

    // True test with equal planes
    Plane p(1, 1, 1, 1), q(1, 1, 1, 1);
    std::cout << equal_plane(p, q) << std::endl;

    // True test with planes whose coefficients are multiple of the other
    Plane pp(1, 1, 1, 1), qq(2, 2, 2, 2);
    std::cout << equal_plane(pp, qq) << std::endl;

    // False test with non equal planes
    Plane p2(1, 2, 1, 1), q2(2, 2, 2, 2);
    std::cout << equal_plane(p2, q2) << std::endl;

    // True test using exact computation
    double x = 1.123456789, y = 1.987654321, e=1e-15;
    Plane ppp(x, x, x, x), qqq(y+e, y+e, y+e, y+e);
    std::cout << equal_plane(ppp, qqq) << std::endl;

    return 0;
}

