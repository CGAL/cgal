#include <cassert>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Voronoi_covariance_3/predicates.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Plane_3 Plane;
typedef CGAL::Voronoi_covariance_3::Plane_dual<K> Plane_dual;

typedef CGAL::Voronoi_covariance_3::Has_on_positive_side_3_dual_point<K>
    Has_on_positive_side_3;

int main () {
    Has_on_positive_side_3 has_on;

    Plane a(1, 1, 0, -1), b(1, 2, 1, -1), c(1, 3, 0, -1); // dual plane : x = 1
    Plane_dual p(a,b,c);
    Plane q(0, 0, 0, -1), qq(2, 0, 0, -1), qqq(1, 1, 0, -1);

    // False
    std::cout << has_on(p, q) << std::endl;
    assert(has_on(p, q) == false);

    // True
    std::cout << has_on(p, qq) << std::endl;
    assert(has_on(p, q) == true);

    // False : on the plane
    std::cout << has_on(p, qqq) << std::endl;
    assert(has_on(p, q) == false);

    return 0;
}

