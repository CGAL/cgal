#include <iostream>
#include <vector>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits_3.h>
#include <CGAL/AABB_segment_primitive_3.h>
#include <cassert>

typedef CGAL::Simple_cartesian<double> K;

typedef K::FT FT;
typedef K::Point_3 Point;
typedef K::Plane_3 Plane;
typedef K::Segment_3 Segment;
typedef K::Triangle_3 Triangle;

typedef std::vector<Segment>::iterator Iterator;
typedef CGAL::AABB_segment_primitive_3<K,Iterator> Primitive;
typedef CGAL::AABB_traits_3<K, Primitive> Traits;
typedef CGAL::AABB_tree<Traits> Tree;

int main()
{
    Point a(1.0, 0.0, 0.0);
    Point b(0.0, 1.0, 0.0);
    Point c(0.0, 0.0, 1.0);
    Point d(0.0, 0.0, 0.0);

    std::vector<Segment> segments;
    segments.push_back(Segment(Point(0, 0, 0),Point(2, 2, 2)));

    Tree tree(segments.begin(),segments.end());
    Plane plane_query(a,b,d);
    Triangle triangle_query(a,b,c);

    // Test calls to all functions
    CGAL::Emptyset_iterator devnull;
    tree.all_intersections(triangle_query, devnull);
    tree.all_intersected_primitives(triangle_query, devnull);
    assert(tree.any_intersected_primitive(triangle_query));
    assert(tree.any_intersection(triangle_query));
    const CGAL::Bbox_3 bbox = tree.bbox();
    assert(bbox == CGAL::Bbox_3(0, 0, 0, 2, 2, 2));
    tree.clear();
    tree.insert(segments.begin(), segments.end());
    tree.build();
    assert(tree.closest_point(Point(-0.1, -0.1, -0.1)) == Point(0, 0, 0));
    assert(tree.closest_point(Point(-0.1, -0.1, -0.1), Point(0, 0, 0)) ==
     Point(0, 0, 0));
    assert(tree.closest_point_and_primitive(Point(-0.1, -0.1, -0.1)).second ==
     segments.begin());
    // Too lazy to call closest_point_and_primitive with a hint. The API is
    // strange. --Laurent Rineau, 2013/01/16
    assert(tree.do_intersect(plane_query) == true);
    assert(tree.do_intersect(triangle_query) == true);
    assert(!tree.empty());
    assert(tree.size() == 1);
    tree.clear();
    assert(tree.size() == 0);
    tree.insert(segments.begin(), segments.end());
    assert(tree.size() == 1);
    assert(tree.number_of_intersected_primitives(plane_query) == 1);
    tree.rebuild(segments.begin(), segments.end());
    assert(tree.size() == 1);
    assert(tree.number_of_intersected_primitives(triangle_query) == 1);
    assert(tree.squared_distance(Point(0, 0, 0)) == 0);

    return EXIT_SUCCESS;
}
