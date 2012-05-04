
// Author(s) : Pierre Alliez

#include <iostream>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/AABB_triangle_primitive.h>
#include <CGAL/Polyhedron_simplex_property_map.h>

typedef CGAL::Simple_cartesian<double> K;
typedef K::FT FT;
typedef K::Point_3 Point;
typedef K::Segment_3 Segment;
typedef CGAL::Polyhedron_3<K> Polyhedron;
typedef CGAL::Triangle_from_facet_property_map<Polyhedron> Triangle_property_map;
typedef CGAL::AABB_triangle_primitive<K,Polyhedron::Facet_handle,Triangle_property_map> Primitive;
typedef CGAL::AABB_traits<K, Primitive> Traits;
typedef CGAL::AABB_tree<Traits> Tree;
typedef Tree::Object_and_primitive_id Object_and_primitive_id;
typedef Tree::Point_and_primitive_id Point_and_primitive_id;

int main()
{
    Point p(1.0, 0.0, 0.0);
    Point q(0.0, 1.0, 0.0);
    Point r(0.0, 0.0, 1.0);
    Point s(0.0, 0.0, 0.0);
    Polyhedron polyhedron;
    polyhedron.make_tetrahedron(p, q, r, s);

    // constructs AABB tree and computes internal KD-tree 
    // data structure to accelerate distance queries
    Tree tree(polyhedron.facets_begin(),polyhedron.facets_end());
    tree.accelerate_distance_queries();

    // query point
    Point query(0.0, 0.0, 3.0);

    // computes squared distance from query
    FT sqd = tree.squared_distance(query);
    std::cout << "squared distance: " << sqd << std::endl;

    // computes closest point
    Point closest = tree.closest_point(query);
    std::cout << "closest point: " << closest << std::endl;

    // computes closest point and primitive id
    Point_and_primitive_id pp = tree.closest_point_and_primitive(query);
    std::cout << "closest point: " << pp.first << std::endl;
    Polyhedron::Face_handle f = pp.second; // closest primitive id
    return EXIT_SUCCESS;
}
