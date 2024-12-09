#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/IO/WKT.h>

#include <iostream>
#include <fstream>
#include <vector>

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
    CGAL::IO::read_WKT(is, points,polylines,polygons);

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
