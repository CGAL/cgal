#include <iostream>
#include <fstream>

#include <CGAL/IO/WKT.h>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#include <boost/foreach.hpp>

#include <vector>

//TODO : If FT is Gmpq, the output of writing will not be doubles.
//typedef CGAL::Simple_cartesian<CGAL::Gmpq> Kernel;

typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;


int main(int argc, char* argv[])
{
  typedef CGAL::Point_2<Kernel> Point;
  typedef CGAL::Geometry_container<std::vector<Point>, boost::geometry::multi_point_tag>  MultiPoint;
  
  std::ifstream is((argc>1)?argv[1]:"data/multipoint.wkt");
  MultiPoint mp;
  CGAL::read_multipoint_WKT(is, mp);
  BOOST_FOREACH(const Point& p, mp)
  {
    std::cout<<p<<std::endl;
  }
}
