// Author(s) : Camille Wormser, Pierre Alliez

#include <iostream>
#include <list>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>

typedef CGAL::Simple_cartesian<double> K;
typedef K::Point_3 Point;
typedef K::Plane_3 Plane;
typedef K::Vector_3 Vector;
typedef K::Segment_3 Segment;
typedef CGAL::Polyhedron_3<K> Polyhedron;
typedef CGAL::AABB_face_graph_triangle_primitive<Polyhedron> Primitive;
typedef CGAL::AABB_traits<K, Primitive> Traits;
typedef CGAL::AABB_tree<Traits> Tree;
typedef boost::optional< Tree::Intersection_and_primitive_id<Segment>::Type > Segment_intersection;
typedef boost::optional< Tree::Intersection_and_primitive_id<Plane>::Type > Plane_intersection;
typedef Tree::Primitive_id Primitive_id;

int main()
{
    Point p(1.0, 0.0, 0.0);
    Point q(0.0, 1.0, 0.0);
    Point r(0.0, 0.0, 1.0);
    Point s(0.0, 0.0, 0.0);
    Polyhedron polyhedron;
    polyhedron.make_tetrahedron(p, q, r, s);

    // constructs AABB tree
    Tree tree(faces(polyhedron).first, faces(polyhedron).second, polyhedron);

    // constructs segment query
    Point a(-0.2, 0.2, -0.2);
    Point b(1.3, 0.2, 1.3);
    Segment segment_query(a,b);

    // tests intersections with segment query
    if(tree.do_intersect(segment_query))
        std::cout << "intersection(s)" << std::endl;
    else
        std::cout << "no intersection" << std::endl;

    // computes #intersections with segment query
    std::cout << tree.number_of_intersected_primitives(segment_query)
        << " intersection(s)" << std::endl;

    // computes first encountered intersection with segment query
    // (generally a point)
    Segment_intersection intersection =
        tree.any_intersection(segment_query);
    if(intersection)
    {
        // gets intersection object
      if(boost::get<Point>(&(intersection->first)))
        std::cout << "intersection object is a point" << std::endl;
    }

    // computes all intersections with segment query (as pairs object - primitive_id)
    std::list<Segment_intersection> intersections;
    tree.all_intersections(segment_query, std::back_inserter(intersections));

    // computes all intersected primitives with segment query as primitive ids
    std::list<Primitive_id> primitives;
    tree.all_intersected_primitives(segment_query, std::back_inserter(primitives));

    // constructs plane query
    Vector vec(0.0,0.0,1.0);
    Plane plane_query(a,vec);

    // computes first encountered intersection with plane query
    // (generally a segment)
    Plane_intersection plane_intersection = tree.any_intersection(plane_query);
    if(plane_intersection)
    {
      
      if(boost::get<Segment>(&(plane_intersection->first)))
            std::cout << "intersection object is a segment" << std::endl;
    }

    return EXIT_SUCCESS;
}
