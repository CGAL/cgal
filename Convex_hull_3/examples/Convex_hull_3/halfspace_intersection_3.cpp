#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Convex_hull_3/dual/halfspace_intersection_3.h>
#include <CGAL/point_generators_3.h>

#include <list>

typedef CGAL::Exact_predicates_inexact_constructions_kernel   K;
typedef K::Plane_3                                            Plane;
typedef K::Point_3                                            Point;
typedef CGAL::Polyhedron_3<K>                                 Polyhedron_3;

// compute the tangent plane of a point
template <typename K>
typename K::Plane_3 tangent_plane (typename K::Point_3 const& p) {
    typename K::Vector_3 v(p.x(), p.y(), p.z());
    v = v / sqrt(v.squared_length());
    typename K::Plane_3 plane(v.x(), v.y(), v.z(), -(p - CGAL::ORIGIN) * v);

    return plane;
}

int main (void) {
    // number of generated planes
    int N = 200;

    // generates random planes on a sphere
    std::list<Plane> planes;
    CGAL::Random_points_on_sphere_3<Point> g;
    for (int i = 0; i < N; i++) {
        planes.push_back(tangent_plane<K>(*g++));
    }

    // define polyhedron to hold the intersection
    Polyhedron_3 P;

    // compute the intersection
    // if no point inside the intersection is provided, one
    // will be automatically found using linear programming
    CGAL::halfspace_intersection_3(planes.begin(),
                                   planes.end(),
                                   P,
                                   boost::make_optional(Point(0, 0, 0)) );

    return 0;
}

