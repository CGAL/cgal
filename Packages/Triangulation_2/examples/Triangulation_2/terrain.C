// file : example/Triangulation_2/terrain.C
#include <CGAL/Homogeneous.h>
#include <fstream>
#include <CGAL/Gmpz.h>
#include <CGAL/Triangulation_euclidean_traits_xy_3.h>
#include <CGAL/Delaunay_triangulation_2.h>

typedef CGAL::Homogeneous<CGAL::Gmpz>  Rp;
typedef CGAL::Triangulation_euclidean_traits_xy_3<Rp>  Gt;
typedef CGAL::Delaunay_triangulation_2<Gt> Delaunay;

typedef Rp::Point_3   Point;

int main()
{
  std::ifstream in("data/terrain.cin");
  std::istream_iterator<Point> begin(in);
  std::istream_iterator<Point> end;

  Delaunay dt;
  dt.insert(begin, end);
  return 0;
}

