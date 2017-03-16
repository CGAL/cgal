#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Periodic_3_regular_triangulation_traits_3.h>
#include <CGAL/Periodic_3_regular_triangulation_3.h>

#include <iostream>

typedef CGAL::Epick K;
typedef K::FT FT;

typedef CGAL::Periodic_3_regular_triangulation_traits_3<K>        Gt;
typedef CGAL::Periodic_3_regular_triangulation_3<Gt>              P3RT3;

typedef Gt::Weighted_point_3                                Weighted_point_3;
typedef Gt::Point_3                                         Point_3;

typedef P3RT3::Vertex_handle                                Vertex_handle;

int main (int, char**)
{
  P3RT3 p3rt3(P3RT3::Iso_cuboid(0,0,0, 1,1,1));

  p3rt3.insert(Weighted_point_3(Point_3(0.1,0.1,0.1), 0.01));
  p3rt3.insert(Weighted_point_3(Point_3(0.9,0.1,0.1), 0.01));
  p3rt3.insert(Weighted_point_3(Point_3(0.1,0.9,0.1), 0.01));
  p3rt3.insert(Weighted_point_3(Point_3(0.1,0.1,0.9), 0.01));
  p3rt3.insert(Weighted_point_3(Point_3(0.4,0.4,0.4), 0.001));

  while (p3rt3.number_of_vertices())
    p3rt3.remove(p3rt3.vertices_begin());

  std::cout << "EXIT SUCCESS" << std::endl;
  return 0;
}
