#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/dual/halfspaces_intersection_with_construction_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>

#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel         K;
typedef K::Plane_3 Plane;
typedef K::Point_3 Point;
typedef CGAL::Polyhedron_3<K> Polyhedron_3;

int main (void) {
    std::list<Plane> planes;
    planes.push_back(Plane(1, 0, 0, -1));
    planes.push_back(Plane(-1, 0, 0, -1));
    planes.push_back(Plane(0, 1, 0, -1));
    planes.push_back(Plane(0, -1, 0, -1));
    planes.push_back(Plane(0, 0, 1, -1));
    planes.push_back(Plane(0, 0, -1, -1));

    Point origin(0, 0, 0);

    Polyhedron_3 P;

    CGAL::halfspaces_intersection_with_construction_3(planes.begin(),
                                                      planes.end(),
                                                      P,
                                                      origin);

    std::ofstream file("halfspaces_intersection_with_construction_3.off");
    file << P;

    return 0;
}

