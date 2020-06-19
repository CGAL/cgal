#include <boost/config.hpp>
#include <boost/version.hpp>

#if BOOST_VERSION >= 105600 && (! defined(BOOST_GCC) || BOOST_GCC >= 40500)
#include <iostream>
#include <fstream>

#include <CGAL/IO/WKT.h>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#include <vector>

//typedef CGAL::Simple_cartesian<CGAL::Gmpq> Kernel;

typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
int main(int argc, char* argv[])
{
  typedef CGAL::Point_2<Kernel> Point;
  typedef std::vector<Point> MultiPoint;

  typedef std::vector<Point> LineString;
  typedef std::vector<LineString> MultiLineString;

  typedef CGAL::Polygon_with_holes_2<Kernel> Polygon;
  typedef std::vector<Polygon> MultiPolygon;

  {
    std::ifstream is((argc>1)?argv[1]:"data/multiple.wkt");
    MultiPoint points;
    MultiLineString polylines;
    MultiPolygon polygons;
    CGAL::read_WKT(is, points,polylines,polygons);

    for(Point p : points)
      std::cout<<p<<std::endl;
    for(LineString ls : polylines)
        for(Point p : ls)
          std::cout<<p<<std::endl;
    for(Polygon p : polygons)
      std::cout<<p<<std::endl;

  }
  return 0;
}
#else
int main()
{
  return 0;
}
#endif
