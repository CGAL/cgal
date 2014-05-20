#include <cassert>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Voronoi_covariance_3/predicates.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Plane_3 Plane;
typedef CGAL::Voronoi_covariance_3::Plane_dual<K> Plane_dual;

typedef CGAL::Voronoi_covariance_3::Less_signed_distance_to_plane_3_dual_point<K>
    Less_signed_distance_to_plane_3;

int main () {
    Less_signed_distance_to_plane_3 less_signed_distance;

    // dual plane : x = 1
    Plane a(1, 1, 0, -1), b(1, 2, 1, -1), c(1, 3, 0, -1);
    Plane_dual p(a,b,c);

    Plane q(2, 0, 0, -1), r(3, 0, 0, -1), rr(0.5, 0, 0, -1);
    // True
    std::cout << less_signed_distance(p, q, r) << std::endl;
    assert(less_signed_distance(p, q, r) == true);

    // False
    std::cout << less_signed_distance(p, q, rr) << std::endl;
    assert(less_signed_distance(p, q, r) == false);

    return 0;
}

