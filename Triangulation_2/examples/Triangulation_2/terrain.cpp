// file: examples/Triangulation_2/terrain.C

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_euclidean_traits_xy_3.h>
#include <CGAL/Delaunay_triangulation_2.h>

#include <fstream>

struct K : CGAL::Exact_predicates_inexact_constructions_kernel {};

typedef CGAL::Triangulation_euclidean_traits_xy_3<K>  Gt;
typedef CGAL::Delaunay_triangulation_2<Gt> Delaunay;

typedef K::Point_3   Point;

int main()
{
  std::ifstream in("data/terrain.cin");
  std::istream_iterator<Point> begin(in);
  std::istream_iterator<Point> end;

  Delaunay dt;
  dt.insert(begin, end);
  std::cout << dt.number_of_vertices() << std::endl;
  return 0;
}

