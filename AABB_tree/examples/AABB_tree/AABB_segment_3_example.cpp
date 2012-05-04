// Author : Pierre Alliez

#include <iostream>
#include <list>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_segment_primitive.h>

typedef CGAL::Simple_cartesian<double> K;

typedef K::FT FT;
typedef K::Point_3 Point;
typedef K::Plane_3 Plane;
typedef K::Segment_3 Segment;
typedef K::Triangle_3 Triangle;

typedef std::list<Segment>::iterator Iterator;
typedef CGAL::AABB_segment_primitive<K,Iterator> Primitive;
typedef CGAL::AABB_traits<K, Primitive> Traits;
typedef CGAL::AABB_tree<Traits> Tree;

int main()
{
    Point a(1.0, 0.0, 0.0);
    Point b(0.0, 1.0, 0.0);
    Point c(0.0, 0.0, 1.0);
    Point d(0.0, 0.0, 0.0);

    std::list<Segment> segments;
    segments.push_back(Segment(a,b));
    segments.push_back(Segment(a,c));
    segments.push_back(Segment(c,d));

    // constructs the AABB tree and the internal search tree for 
    // efficient distance computations.
    Tree tree(segments.begin(),segments.end());
    tree.accelerate_distance_queries();

    // counts #intersections with a plane query
    Plane plane_query(a,b,d);
    std::cout << tree.number_of_intersected_primitives(plane_query)
        << " intersections(s) with plane" << std::endl;

    // counts #intersections with a triangle query
    Triangle triangle_query(a,b,c);
    std::cout << tree.number_of_intersected_primitives(triangle_query)
        << " intersections(s) with triangle" << std::endl;

    // computes the closest point from a point query 
    Point point_query(2.0, 2.0, 2.0);
    Point closest = tree.closest_point(point_query);

    std::cerr << "closest point is: " << closest << std::endl;
    return EXIT_SUCCESS;
}
