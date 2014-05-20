#include <cassert>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Voronoi_covariance_3/predicates.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Plane_3 Plane;

typedef CGAL::Voronoi_covariance_3::Coplanar_3_dual_point<K> Coplanar_3;

int main () {
    Coplanar_3 coplanar;

    // Simple true test with equal planes
    K::Plane_3 p(1, 1, 1, 1), q(1, 1, 1, 1), r(1, 1, 1, 1), s(1, 1, 1, 1);
    std::cout << coplanar(p, q, r, s) << std::endl;
    assert(coplanar(p, q, r, s) == true);

    // Simple true test
    K::Plane_3 p2(1, 1, 1, -1), q2(2, 1, 1, -1), r2(2, 2, 1, -1), s2(1, 2, 1, -1);
    std::cout << coplanar(p2, q2, r2, s2) << std::endl;
    assert(coplanar(p2, q2, r2, s2) == true);

    // Simple false test
    K::Plane_3 p3(1, 1, 1, -1), q3(2, 1, 1, -1), r3(2, 2, 1, -1), s3(1, 1, 2, -1);
    std::cout << coplanar(p3, q3, r3, s3) << std::endl;
    assert(coplanar(p3, q3, r3, s3) == false);

    return 0;
}

