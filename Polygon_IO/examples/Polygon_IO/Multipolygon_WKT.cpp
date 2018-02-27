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
  typedef CGAL::Polygon_with_holes_2<Kernel> Polygon;
  typedef CGAL::Point_2<Kernel> Point;
  std::ifstream is((argc>1)?argv[1]:"data/multipolygon.wkt");
  std::vector<Polygon> multi_poly;
  CGAL::read_multipolygon_WKT(is, multi_poly);
  
  BOOST_FOREACH(Polygon p, multi_poly)
      std::cout<<p<<std::endl;
  
}
