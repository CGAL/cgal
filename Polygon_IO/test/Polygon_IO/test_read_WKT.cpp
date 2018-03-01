#include <iostream>
#include <fstream>

#include <CGAL/IO/WKT.h>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <boost/foreach.hpp>

#include <vector>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;


int main()
{
  typedef CGAL::Point_2<Kernel> Point;
  typedef CGAL::Geometry_container<std::vector<Point>, boost::geometry::linestring_tag> Linestring;
  typedef CGAL::Polygon_with_holes_2<Kernel> Polygon;
  typedef CGAL::Geometry_container<std::vector<Point>, boost::geometry::multi_point_tag> MultiPoint;
  typedef CGAL::Geometry_container<std::vector<Linestring>, boost::geometry::multi_linestring_tag> MultiLinestring;
  typedef CGAL::Geometry_container<std::vector<Polygon>, boost::geometry::multi_polygon_tag> MultiPolygon;
  
  Point p;
  std::ifstream in("data/point.wkt");
  CGAL::read_point_WKT(in, p);
  in.close();
  CGAL_assertion(p == Point(2,3));
  
  in = std::ifstream("data/linestring.wkt");
  Linestring ls;
  CGAL::read_linestring_WKT(in, ls);
  in.close();
  CGAL_assertion(ls.size() == 3);
  
  Polygon poly;
  in = std::ifstream("data/polygon.wkt");
  CGAL::read_polygon_WKT(in, poly);
  in.close();
  CGAL_assertion(poly.outer_boundary().size() == 3);
  
  
  MultiPoint pees;
  in = std::ifstream("data/multipoint.wkt");
  CGAL::read_multipoint_WKT(in, pees);
  in.close();
  CGAL_assertion(pees.size() == 4);
  
  in = std::ifstream("data/multilinestring.wkt");
  MultiLinestring mls;
  CGAL::read_multilinestring_WKT(in, mls);
  in.close();
  CGAL_assertion(mls.size() == 2);
  
  MultiPolygon polies;
  in = std::ifstream("data/multipolygon.wkt");
  CGAL::read_multipolygon_WKT(in, polies);
  in.close();
  CGAL_assertion(polies.size() == 2);
  return 0;
}
