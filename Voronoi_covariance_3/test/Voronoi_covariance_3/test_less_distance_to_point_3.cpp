#include <cassert>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Voronoi_covariance_3/predicates.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Plane_3 Plane;

typedef CGAL::Voronoi_covariance_3::Less_distance_to_point_3_dual_point<K>
    Less_distance_to_point_3;

int main () {
    Less_distance_to_point_3 less_distance_to_point;

    Plane p(1, 0, 0, -1), q(3, 0, 0, -1);
    Plane r(2, 0, 0, -1), rr(4, 0, 0, -1);

    // True
    std::cout << less_distance_to_point(p, q, rr) << std::endl;
    assert(less_distance_to_point(p, q, rr) == true);

    // False
    std::cout << less_distance_to_point(p, q, r) << std::endl;
    assert(less_distance_to_point(p, q, r) == false);

    // Tests if r is one of the two first points
    // False
    std::cout << less_distance_to_point(p, q, p) << std::endl;
    assert(less_distance_to_point(p, q, p) == false);
    std::cout << less_distance_to_point(p, q, q) << std::endl;
    assert(less_distance_to_point(p, q, q) == false);

    return 0;
}

