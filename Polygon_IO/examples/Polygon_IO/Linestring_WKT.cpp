#include <iostream>
#include <fstream>
#include <CGAL/IO/WKT.h>
#include <boost/foreach.hpp> //must be included before WKT for some reason
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <vector>

//TODO : If FT is Gmpq, the output of writing will not be doubles.
//typedef CGAL::Simple_cartesian<CGAL::Gmpq> Kernel;

typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;


int main(int argc, char* argv[])
{
  typedef CGAL::Point_2<Kernel> Point;
  typedef CGAL::Geometry_container<std::vector<Point>,boost::geometry::linestring_tag> LineString;
  typedef CGAL::Geometry_container<std::vector<LineString>, boost::geometry::multi_linestring_tag> MultiLineString;
  LineString ls;
  {
    std::ifstream is((argc>1)?argv[1]:"data/linestring.wkt");
    //std::vector<Point> ls;
    CGAL::read_linestring_WKT(is, ls);
    is.close();
  }
  BOOST_FOREACH(Point p, ls)
      std::cout<<p<<std::endl;
  ls.clear();
  MultiLineString mls;
  {
    std::ifstream is((argc>2)?argv[2]:"data/multilinestring.wkt");
    CGAL::read_multilinestring_WKT(is, mls);
    is.close();
  }
  BOOST_FOREACH(LineString l, mls)
  {
    BOOST_FOREACH(const Point& p, l)
        std::cout<<p<<std::endl;
  }
  
}
