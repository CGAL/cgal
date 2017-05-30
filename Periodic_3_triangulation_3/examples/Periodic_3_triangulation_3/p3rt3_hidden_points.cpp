#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Periodic_3_regular_triangulation_traits_3.h>
#include <CGAL/Periodic_3_regular_triangulation_3.h>

#include <iostream>

typedef CGAL::Epick                                              K;
typedef CGAL::Periodic_3_regular_triangulation_traits_3<K>       Gt;
typedef CGAL::Periodic_3_regular_triangulation_3<Gt>             P3RT3;

typedef Gt::Iso_cuboid_3                                     Iso_cuboid;
typedef Gt::Point_3                                          Point_3;
typedef Gt::Weighted_point_3                                 Weighted_point_3;

typedef P3RT3::Vertex_handle                                 Vertex_handle;

int main ()
{
  P3RT3 p3rt3(P3RT3::Iso_cuboid(0,0,0, 1,1,1));

  p3rt3.insert(Weighted_point_3(Point_3(0.1,0.1,0.1), 0.01));
  p3rt3.insert(Weighted_point_3(Point_3(0.9,0.1,0.1), 0.01));
  p3rt3.insert(Weighted_point_3(Point_3(0.1,0.91,0.1), 0.01));
  p3rt3.insert(Weighted_point_3(Point_3(0.1,0.1,0.89), 0.01));
  p3rt3.insert(Weighted_point_3(Point_3(0.101,0.101,0.101), 0.001));
  assert(p3rt3.is_valid());

  std::cout << "Number of vertices :      " << p3rt3.number_of_vertices() << std::endl;
  std::cout << "Number of hidden points : " << p3rt3.number_of_hidden_points() << std::endl;

  std::cout << "Removing the first point..." << std::endl;
  p3rt3.remove(p3rt3.vertices_begin());
  // The first point is removed and the hidden point is revealed.

  std::cout << "Number of vertices :      " << p3rt3.number_of_vertices() << std::endl;
  std::cout << "Number of hidden points : " << p3rt3.number_of_hidden_points() << std::endl;

  std::cout << "EXIT SUCCESS" << std::endl;
  return 0;
}
