#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/dual/halfspaces_intersection.h>

#include <list>

typedef CGAL::Exact_predicates_inexact_constructions_kernel   K;
typedef K::Plane_3                                            Plane;
typedef CGAL::Polyhedron_3<K>                                 Polyhedron_3;

int main (void) {
    // planes corresponding to the faces of a cube
    std::list<Plane> planes;
    planes.push_back(Plane(1, 0, 0, -1));
    planes.push_back(Plane(-1, 0, 0, -1));
    planes.push_back(Plane(0, 1, 0, -1));
    planes.push_back(Plane(0, -1, 0, -1));
    planes.push_back(Plane(0, 0, 1, -1));
    planes.push_back(Plane(0, 0, -1, -1));

    // define polyhedron to hold the intersection
    Polyhedron_3 P;

    // compute the intersection
    CGAL::halfspaces_intersection(planes.begin(),
                                  planes.end(),
                                  P);

    return 0;
}

