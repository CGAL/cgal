#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Projection_traits_xy_3.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>

#include <string>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Projection_traits_xy_3<K>  Gt;
typedef CGAL::Triangulation_vertex_base_with_info_2<std::string, Gt> Vb;
typedef CGAL::Triangulation_data_structure_2<Vb> Tds;
typedef CGAL::Delaunay_triangulation_2<Gt,Tds> Delaunay;

typedef K::Point_3   Point;

int main()
{
  Delaunay dt;
  Delaunay::Vertex_handle vh;

  vh  = dt.insert(Point(0,0,0));
  vh->info() = "Paris";
  vh = dt.insert(Point(1,0,0.1));
  vh->info() = "London";
  vh = dt.insert(Point(0,1,0.2));
  vh->info() = "New York";

  return 0;
}
