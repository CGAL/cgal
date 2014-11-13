#include <iostream>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>


typedef CGAL::Simple_cartesian<double> K;
typedef K::FT FT;
typedef K::Point_3 Point;
typedef K::Segment_3 Segment;
typedef CGAL::Polyhedron_3<K> Polyhedron;
typedef CGAL::AABB_face_graph_triangle_primitive<Polyhedron, CGAL::Default, CGAL::Tag_false> Primitive;
typedef CGAL::AABB_traits<K, Primitive> Traits;
typedef CGAL::AABB_tree<Traits> Tree;
typedef Tree::Point_and_primitive_id Point_and_primitive_id;

int main()
{
    Point p(1.0, 0.0, 0.0);
    Point q(0.0, 1.0, 0.0);
    Point r(0.0, 0.0, 1.0);
    Point s(0.0, 0.0, 0.0);
    Polyhedron polyhedron1;
    polyhedron1.make_tetrahedron(p, q, r, s);


    Point p2(11.0, 0.0, 0.0);
    Point q2(10.0, 1.0, 0.0);
    Point r2(10.0, 0.0, 1.0);
    Point s2(10.0, 0.0, 0.0);
    Polyhedron polyhedron2;
    polyhedron2.make_tetrahedron(p2, q2, r2, s2);
    // constructs AABB tree and computes internal KD-tree
    // data structure to accelerate distance queries
    Tree tree(faces(polyhedron1).first, faces(polyhedron1).second, polyhedron1);

    tree.accelerate_distance_queries();

    tree.insert(faces(polyhedron2).first, faces(polyhedron2).second, polyhedron2);

    // query point
    Point query(0.0, 0.0, 3.0);

    // computes squared distance from query
    FT sqd = tree.squared_distance(query);
    std::cout << "squared distance: " << sqd << std::endl;

    return EXIT_SUCCESS;
}
