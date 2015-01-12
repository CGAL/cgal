// Author(s) : Pierre Alliez

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
typedef CGAL::AABB_face_graph_triangle_primitive<Polyhedron> Primitive;
typedef CGAL::AABB_traits<K, Primitive> Traits;
typedef CGAL::AABB_tree<Traits> Tree;

template <class Kernel, class FaceGraph>
void run(const FaceGraph& graph){
  typename Kernel::Point_3 p(0.2, 0.2, 0.2);
  typename Kernel::Point_3 q(1.0, 1.0, 1.0);

  // constructs the AABB tree and the internal search tree for
  // efficient distance queries.
  Tree tree( faces(graph).first, faces(graph).second, graph);
  tree.accelerate_distance_queries();

  // counts #intersections with a triangle query
  Segment segment_query(p,q);
  std::cout << tree.number_of_intersected_primitives(segment_query)
      << " intersections(s) with triangle" << std::endl;

  assert( tree.number_of_intersected_primitives(segment_query)== 1 );

  // computes the closest point from a query point
  typename Kernel::Point_3 point_query(2.0, 2.0, 2.0);
  typename Kernel::Point_3 closest = tree.closest_point(point_query);

  std::cerr << "closest point is: " << closest << std::endl;
}

int main()
{
  Point p(1.0, 0.0, 0.0);
  Point q(0.0, 1.0, 0.0);
  Point r(0.0, 0.0, 1.0);
  Point s(0.0, 0.0, 0.0);
  Polyhedron polyhedron;
  polyhedron.make_tetrahedron(p, q, r, s);

  run<K>(polyhedron);
  return EXIT_SUCCESS;
}
