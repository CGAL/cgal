#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/IO/WKT.h>

#include <iostream>
#include <fstream>
#include <vector>

typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;

int main(int argc, char* argv[])
{
  typedef CGAL::Point_2<Kernel> Point;
  typedef std::vector<Point>  MultiPoint;

  std::ifstream is((argc>1)?argv[1]:"data/multipoint.wkt");
  MultiPoint mp;
  CGAL::IO::read_multi_point_WKT(is, mp);
  for(const Point& p : mp)
  {
    std::cout<<p<<std::endl;
  }
  is.close();
  return 0;
}
