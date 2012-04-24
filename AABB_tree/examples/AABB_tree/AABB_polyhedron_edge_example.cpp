// Author(s) : Pierre Alliez

#include <iostream>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/AABB_polyhedron_segment_primitive.h>

typedef CGAL::Simple_cartesian<double> K;
typedef K::FT FT;
typedef K::Point_3 Point;
typedef K::Triangle_3 Triangle;
typedef CGAL::Polyhedron_3<K> Polyhedron;
typedef CGAL::AABB_polyhedron_segment_primitive<K,Polyhedron> Primitive;
typedef CGAL::AABB_traits<K, Primitive> Traits;
typedef CGAL::AABB_tree<Traits> Tree;

int main()
{
    Point p(1.0, 0.0, 0.0);
    Point q(0.0, 1.0, 0.0);
    Point r(0.0, 0.0, 1.0);
    Point s(0.0, 0.0, 0.0);
    Polyhedron polyhedron;
    polyhedron.make_tetrahedron(p, q, r, s);

    // constructs the AABB tree and the internal search tree for 
    // efficient distance queries.
    Tree tree(polyhedron.edges_begin(),polyhedron.edges_end());
    tree.accelerate_distance_queries();

    // counts #intersections with a triangle query
    Triangle triangle_query(p,q,r);
    std::cout << tree.number_of_intersected_primitives(triangle_query)
        << " intersections(s) with triangle" << std::endl;

    // computes the closest point from a query point
    Point point_query(2.0, 2.0, 2.0);
    Point closest = tree.closest_point(point_query);

    std::cerr << "closest point is: " << closest << std::endl;
    return EXIT_SUCCESS;
}
