#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/dual/halfspaces_intersection_3.h>

#include "include/off.h"

typedef CGAL::Exact_predicates_inexact_constructions_kernel         K;
typedef K::Plane_3 Plane;
typedef K::Point_3 Point;
typedef K::RT RT;
typedef CGAL::Polyhedron_3<K> Polyhedron_3;

// Translate a plane
Plane translate_plane (Plane const& p, Point const& p0) {
    RT newA = p.a();
    RT newB = p.b();
    RT newC = p.c();
    RT newD = p.d() + p.a() * p0.x() + p.b() * p0.y() + p.c() * p0.z();

    Plane newP(newA, newB, newC, newD);
    return newP;
}

#include <CGAL/point_generators_3.h>

// Tangent plane
template <typename K>
typename K::Plane_3 tangent_plane (typename K::Point_3 const& p) {
    typename K::Vector_3 v(p.x(), p.y(), p.z());
    v = v / sqrt(v.squared_length());
    typename K::Plane_3 plane(v.x(), v.y(), v.z(), -(p - CGAL::ORIGIN) * v);

    return plane;
}

int main (int argc, char *argv[]) {
    // Translated cuboid
    Point o(0.1, 1.5, -0.2);
    std::list<Plane> planes;
    planes.push_back(Plane(1, 0, 0, -1));
    planes.push_back(Plane(-1, 0, 0, -1));

    planes.push_back(Plane(0, 1, 0, -2));
    planes.push_back(Plane(0, 1, 0, -1));

    planes.push_back(Plane(0, 0, 1, -1));
    planes.push_back(Plane(0, 0, -1, -1));

    // Random points on a translated sphere
    std::list<Plane> sphere_planes;
    Point os(0, 0, 2), oos(0, 0, -2);
    int N;
    if (argc > 1) {
        N = atoi(argv[1]);
    } else {
        N = 200;
    }

    CGAL::Random_points_on_sphere_3<Point> g;
    for (int i = 0; i < N; i++) {
        Plane p = tangent_plane<K>(*g++);
        Plane pp = translate_plane(p, oos);
        sphere_planes.push_back(pp);
    }

    Polyhedron_3 P, PP;

    CGAL::halfspaces_intersection_3(planes.begin(), planes.end(), P, o);
    convertToOFF<K, Polyhedron_3>("translated_cube.off", P);

    CGAL::halfspaces_intersection_3(sphere_planes.begin(), sphere_planes.end(), PP, os);
    convertToOFF<K, Polyhedron_3>("translated_sphere.off", PP);

    return 0;
}

