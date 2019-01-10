#include <iostream>
#include <fstream>
#if BOOST_VERSION >= 105600 && (! defined(BOOST_GCC) || BOOST_GCC >= 40500)
#include <CGAL/IO/WKT.h>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <boost/foreach.hpp>

#include <vector>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;


int main()
{
  typedef CGAL::Point_2<Kernel> Point;
  typedef std::vector<Point> Linestring;
  typedef CGAL::Polygon_with_holes_2<Kernel> Polygon;
  typedef std::vector<Point> MultiPoint;
  typedef std::vector<Linestring> MultiLinestring;
  typedef std::vector<Polygon> MultiPolygon;
  
  Point p;
  {
    std::ifstream in("data/point.wkt");
    CGAL::read_point_WKT(in, p);
    in.close();
    CGAL_assertion(p == Point(2,3));
  }
  {
    std::ifstream in("data/linestring.wkt");
    Linestring ls;
    CGAL::read_linestring_WKT(in, ls);
    in.close();
    CGAL_assertion(ls.size() == 3);
  }
  {
    Polygon poly;
    std::ifstream in("data/polygon.wkt");
    CGAL::read_polygon_WKT(in, poly);
    in.close();
    CGAL_assertion(poly.outer_boundary().size() == 3);
  }
  {
    MultiPoint pees;
    std::ifstream in("data/multipoint.wkt");
    CGAL::read_multi_point_WKT(in, pees);
    in.close();
    CGAL_assertion(pees.size() == 4);
  }
  {
    std::ifstream in("data/multilinestring.wkt");
    MultiLinestring mls;
    CGAL::read_multi_linestring_WKT(in, mls);
    in.close();
    CGAL_assertion(mls.size() == 2);
  }
  {
    MultiPolygon polies;
    std::ifstream in("data/multipolygon.wkt");
    CGAL::read_multi_polygon_WKT(in, polies);
    in.close();
    CGAL_assertion(polies.size() == 2);
  }
  std::cout<<"WKT reading test passed."<<std::endl;
  return 0;
}
#else
int main()
{
  return 0;
}
#endif
