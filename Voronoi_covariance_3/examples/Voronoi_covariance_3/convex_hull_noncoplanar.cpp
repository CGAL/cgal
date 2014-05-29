#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Polyhedron_3.h>
#include <CGAL/Voronoi_covariance_3/Convex_hull_traits_dual_3.h>

// TODO: primal polyhedron

typedef CGAL::Exact_predicates_inexact_constructions_kernel         K;
typedef K::Plane_3 Plane;
typedef K::Point_3 Point;

typedef CGAL::Voronoi_covariance_3::Convex_hull_traits_dual_3<K>    Hull_traits_dual_3;
typedef CGAL::Polyhedron_3<Hull_traits_dual_3>                      Polyhedron_dual_3;

// Specialization for dual traits
namespace CGAL {
    namespace internal {
        namespace Convex_hull_3 {
            template <class InputIterator, class Plane_3, class Polyhedron_3>
                void coplanar_3_hull(InputIterator first, InputIterator beyond,
                                     Plane_3 plane, Polyhedron_3& P, const Hull_traits_dual_3& traits) {
                    // do nothing
                }
        }
    }
}

#include <CGAL/convex_hull_3.h>
#include "include/to_dual.h"

#include <CGAL/point_generators_3.h>

template <typename K>
typename K::Plane_3 tangent_plane (typename K::Point_3 const& p) {
    typename K::Vector_3 v(p.x(), p.y(), p.z());
    typename K::Plane_3 plane(v.x(), v.y(), v.z(), -(p - CGAL::ORIGIN) * v);
    v = v / sqrt(v.squared_length());

    return plane;
}

int main (void) {
    // define polyhedron to hold convex hull
    Polyhedron_dual_3 dual_poly_sphere;
    Polyhedron_dual_3 dual_poly_cube;

    // traits
    Hull_traits_dual_3 dual_traits;

    // Cube
    // IMPORTANT: d <= 0
    std::list<Plane> planes;
    planes.push_back(Plane(1, 0, 0, -1));
    planes.push_back(Plane(-1, 0, 0, -1));
    planes.push_back(Plane(0, 1, 0, -1));
    planes.push_back(Plane(0, -1, 0, -1));
    planes.push_back(Plane(0, 0, 1, -1));
    planes.push_back(Plane(0, 0, -1, -1));

    // Random points on a sphere
    std::vector<Plane> sphere_planes;
    int N = 200;

    CGAL::Random_points_on_sphere_3<Point> g;
    for (int i = 0; i < N; i++) {
        sphere_planes.push_back(tangent_plane<K>(*g++));
    }

    // Compute dual convex hull
    CGAL::convex_hull_3(planes.begin(), planes.end(), dual_poly_cube, dual_traits);
    CGAL::convex_hull_3(sphere_planes.begin(), sphere_planes.end(), dual_poly_sphere, dual_traits);

    // Print the polyhedron in an OFF file
    convert_dual_OFF<K>("dual_cube_convex_hull.off", dual_poly_cube);
    convert_dual_OFF<K>("dual_sphere_convex_hull.off", dual_poly_sphere);

    return 0;
}

