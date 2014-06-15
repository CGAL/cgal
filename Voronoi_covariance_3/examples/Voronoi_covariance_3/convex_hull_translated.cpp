#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Voronoi_covariance_3/halfspaces_intersection.h>

#include "include/off.h"

typedef CGAL::Exact_predicates_inexact_constructions_kernel         K;
typedef K::Plane_3 Plane;
typedef K::Point_3 Point;
typedef CGAL::Polyhedron_3<K> Polyhedron_3;
typedef K::RT RT;

Plane translate_plane (Plane const& p, Point const& p0) {
    RT newA = p.a();
    RT newB = p.b();
    RT newC = p.c();
    RT newD = p.d() + p.a() * p0.x() + p.b() * p0.y() + p.c() * p0.z();

    Plane newP(newA, newB, newC, newD);

    return newP;
}

int main (void) {
    std::list<Plane> planes;
    planes.push_back(Plane(1, 0, 0, -1));
    planes.push_back(Plane(-1, 0, 0, -1));

    planes.push_back(Plane(0, 1, 0, -2));
    planes.push_back(Plane(0, -1, 0, -1));

    planes.push_back(Plane(0, 0, 1, -2));
    planes.push_back(Plane(0, 0, -1, -1));

    Point o(0, 1, 0);

    std::list<Plane> translated_planes;
    for (std::list<Plane>::const_iterator it = planes.begin();
         it != planes.end();
         it++) {
        translated_planes.push_back(translate_plane(*it, o));
    }

    Polyhedron_3 P;

    /* CGAL::Voronoi_covariance_3::halfspaces_intersection(translated_planes.begin(), translated_planes.end(), P, K(), o); */
    CGAL::Voronoi_covariance_3::halfspaces_intersection(planes.begin(), planes.end(), P, K(), o);

    convertToOFF<K, Polyhedron_3>("translated_cube.off", P);

    return 0;
}

