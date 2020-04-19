#include <iostream>
#include <fstream>
#include <boost/config.hpp>
#include <boost/version.hpp>

#if BOOST_VERSION >= 105600 && (! defined(BOOST_GCC) || BOOST_GCC >= 40500)
#include <CGAL/IO/WKT.h>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#include <vector>

//typedef CGAL::Simple_cartesian<CGAL::Gmpq> Kernel;

typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;


int main(int argc, char* argv[])
{
  typedef CGAL::Point_2<Kernel> Point;
  typedef std::vector<Point>  MultiPoint;

  std::ifstream is((argc>1)?argv[1]:"data/multipoint.wkt");
  MultiPoint mp;
  CGAL::read_multi_point_WKT(is, mp);
  for(const Point& p : mp)
  {
    std::cout<<p<<std::endl;
  }
  is.close();
  return 0;
}
#else
int main()
{
  return 0;
}
#endif
