#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Projection_traits_xy_3.h>
#include <CGAL/Delaunay_triangulation_2.h>

#include <fstream>


typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Projection_traits_xy_3<K>  Gt;
typedef CGAL::Delaunay_triangulation_2<Gt> Delaunay;

typedef K::Point_3   Point;

int main()
{
  std::ifstream in("data/terrain.cin");
  std::istream_iterator<Point> begin(in);
  std::istream_iterator<Point> end;

  Delaunay dt(begin, end);
  std::cout << dt.number_of_vertices() << std::endl;
  return 0;
}
