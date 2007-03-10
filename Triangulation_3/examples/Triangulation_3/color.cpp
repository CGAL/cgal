#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <CGAL/IO/Color.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

typedef CGAL::Triangulation_vertex_base_with_info_3<CGAL::Color, K> Vb;
typedef CGAL::Triangulation_data_structure_3<Vb>                    Tds;
typedef CGAL::Delaunay_triangulation_3<K, Tds>                      Delaunay;

typedef Delaunay::Point   Point;

int main()
{
  Delaunay T;

  T.insert(Point(0,0,0));
  T.insert(Point(1,0,0));
  T.insert(Point(0,1,0));
  T.insert(Point(0,0,1));
  T.insert(Point(2,2,2));
  T.insert(Point(-1,0,1));

  // Set the color of finite vertices of degree 6 to red.
  Delaunay::Finite_vertices_iterator vit;
  for (vit = T.finite_vertices_begin(); vit != T.finite_vertices_end(); ++vit)
      if (T.degree(vit) == 6)
          vit->info() = CGAL::RED;

  return 0;
}
