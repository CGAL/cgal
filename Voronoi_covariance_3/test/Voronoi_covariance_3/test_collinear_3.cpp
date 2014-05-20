#include <cassert>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Voronoi_covariance_3/predicates.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Plane_3 Plane;

typedef CGAL::Voronoi_covariance_3::Collinear_3_dual_point<K> Collinear_3_dual;

int main () {
    Collinear_3_dual collinear;

    // Simple true test with equal planes
    Plane p(1, 1, 1, 1), q(1, 1, 1, 1), r(1, 1, 1, 1);
    std::cout << collinear(p, q, r) << std::endl;
    assert(collinear(p, q, r) == true);

    // Simple true test with non equal planes
    Plane p2(1, 1, 1, -1), q2(2, 2, 2, -1), r2(3, 3, 3, -1);
    std::cout << collinear(p2, q2, r2) << std::endl;
    assert(collinear(p2, q2, r2) == true);

    // Simple false test
    Plane p3(1, 1, 1, -1), q3(2, 2, 2, -1), r3(3, 3, 2, -1);
    std::cout << collinear(p3, q3, r3) << std::endl;
    assert(collinear(p3, q3, r3) == false);

    // Almost but not collinear
    double almost1 = 0.9999999999;
    Plane p4(0, 0, 0, -1), q4(1, 1, 1, -1);
    Plane r4(almost1, almost1, almost1+1e-15, -1);
    std::cout << collinear(p4, q4, r4) << std::endl;
    assert(collinear(p4, q4, r4) == false);

    return 0;
}

