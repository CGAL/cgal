#include <iostream>
#include <vector>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_segment_primitive.h>
#include <cassert>

typedef CGAL::Simple_cartesian<double> K;

typedef K::FT FT;
typedef K::Point_3 Point;
typedef K::Plane_3 Plane;
typedef K::Segment_3 Segment;
typedef K::Triangle_3 Triangle;

typedef std::vector<Segment>::iterator Iterator;
typedef CGAL::AABB_segment_primitive<K,Iterator> Primitive;
typedef CGAL::AABB_traits<K, Primitive> Traits;
typedef CGAL::AABB_tree<Traits> Tree;

int main()
{
    Point a(1.0, 0.0, 0.0);
    Point b(0.0, 1.0, 0.0);
    Point c(0.0, 0.0, 1.0);
    Point d(0.0, 0.0, 0.0);

    Tree tree;
    Plane plane_query(a,b,d);
    Triangle triangle_query(a,b,c);

    // Test calls to all functions but those who have `!empty()` as
    // precondition.
    CGAL::Emptyset_iterator devnull;
    tree.all_intersections(triangle_query, devnull);
    tree.all_intersected_primitives(triangle_query, devnull);
    assert(!tree.any_intersected_primitive(triangle_query));
    assert(!tree.any_intersection(triangle_query));
    //Cannot call tree.bbox();
    tree.build();
    tree.clear();
    //Cannot call tree.closest_*(...)
    assert(tree.do_intersect(plane_query) == false);
    assert(tree.do_intersect(triangle_query) == false);
    assert(tree.empty());
    //Do not call tree.insert(...)
    assert(tree.number_of_intersected_primitives(plane_query) == 0);
    assert(tree.number_of_intersected_primitives(triangle_query) == 0);
    // Cannot call tree.rebuild(..)
    assert(tree.size() == 0);
    // Cannot call tree.squared_distance(..)

    return EXIT_SUCCESS;
}
