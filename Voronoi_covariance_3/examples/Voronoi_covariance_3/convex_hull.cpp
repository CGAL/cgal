#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Polyhedron_3.h>
#include <CGAL/Convex_hull_3/dual/Convex_hull_traits_dual_3.h>
#include <CGAL/Convex_hull_3/dual/halfspace_intersection_3.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel         K;
typedef K::Plane_3 Plane;
typedef K::Point_3 Point;
typedef CGAL::Polyhedron_3<K> Polyhedron_3;

typedef CGAL::Convex_hull_3::Convex_hull_traits_dual_3<K>    Hull_traits_dual_3;
typedef CGAL::Polyhedron_3<Hull_traits_dual_3>                      Polyhedron_dual_3;

#include <CGAL/convex_hull_3.h>
#include "include/to_dual.h"
#include "include/off.h"

#include <CGAL/point_generators_3.h>

template <typename K>
typename K::Plane_3 tangent_plane (typename K::Point_3 const& p) {
    typename K::Vector_3 v(p.x(), p.y(), p.z());
    v = v / sqrt(v.squared_length());
    typename K::Plane_3 plane(v.x(), v.y(), v.z(), -(p - CGAL::ORIGIN) * v);

    return plane;
}

int main (int argc, char *argv[]) {
    // define dual polyhedrons to hold convex hulls
    Polyhedron_dual_3 dual_poly_sphere;
    Polyhedron_dual_3 dual_poly_cube;

    // primal polyhedrons
    Polyhedron_3 poly_cube;
    Polyhedron_3 poly_sphere;

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
    std::list<Plane> sphere_planes;
    int N;
    if (argc > 1) {
        N = atoi(argv[1]);
    } else {
        N = 200;
    }

    CGAL::Random_points_on_sphere_3<Point> g;
    for (int i = 0; i < N; i++) {
        sphere_planes.push_back(tangent_plane<K>(*g++));
    }

    // Compute dual convex hulls
    CGAL::convex_hull_3(planes.begin(), planes.end(), dual_poly_cube, dual_traits);
    CGAL::convex_hull_3(sphere_planes.begin(), sphere_planes.end(), dual_poly_sphere, dual_traits);

    // Print dual polyhedrons in an OFF file
    convert_dual_OFF<K>("dual_cube_convex_hull.off", dual_poly_cube);
    convert_dual_OFF<K>("dual_sphere_convex_hull.off", dual_poly_sphere);

    // Compute associated primal polyhedrons
    typedef CGAL::Convex_hull_3::internal::Build_primal_polyhedron<K, Polyhedron_dual_3, Polyhedron_3>
        Build_primal_polyhedron;
    Build_primal_polyhedron bpp_cube(dual_poly_cube);
    Build_primal_polyhedron bpp_sphere(dual_poly_sphere);

    poly_cube.delegate(bpp_cube);
    poly_sphere.delegate(bpp_sphere);

    // Print primal polyhedrons into an OFF file
    convertToOFF<K, Polyhedron_3>("primal_cube.off", poly_cube);
    convertToOFF<K, Polyhedron_3>("primal_sphere.off", poly_sphere);

    return 0;
}

