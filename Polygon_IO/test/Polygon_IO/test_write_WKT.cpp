#include <iostream>
#include <fstream>

#include <CGAL/IO/WKT.h>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <boost/foreach.hpp>

#include <vector>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Point_2<Kernel> Point;
typedef CGAL::Geometry_container<std::vector<Point>, boost::geometry::linestring_tag> Linestring;
typedef CGAL::Polygon_with_holes_2<Kernel> Polygon;
typedef CGAL::Geometry_container<std::vector<Point>, boost::geometry::multi_point_tag> MultiPoint;
typedef CGAL::Geometry_container<std::vector<Linestring>, boost::geometry::multi_linestring_tag> MultiLinestring;
typedef CGAL::Geometry_container<std::vector<Polygon>, boost::geometry::multi_polygon_tag> MultiPolygon;

double fRand(double fMin, double fMax)
{
  double f = (double)rand() / RAND_MAX;
  return fMin + f * (fMax - fMin);
}

Point generate_point(double xmin, double xmax,
                     double ymin, double ymax)
{
  double x= fRand(xmin,xmax),
      y = fRand(ymin,ymax);
  return Point(x,y);
}

Linestring generate_linestring()
{
  Linestring ls;
  ls.push_back(generate_point(0,15,0,15));
  ls.push_back(generate_point(0,15,0,15));
  return ls;
}
Polygon generate_polygon()
{
  Point bl,br, t;
  bl = generate_point(-10,-5, -10, -5);
  br = generate_point(5,10,-10,-5);
  t = generate_point(-4.99,4.99,5,10);
  br = Point(br.x(), bl.y());
  
  double xmax(br.x()),
      ymax(t.y()),
      xmin(bl.x()),
      xt(t.x()),
      ymin(bl.y()),
      ymid((ymax+ymin)/4.0);
  Polygon::Polygon_2 border;
  border.push_back(bl);
  border.push_back(t);
  border.push_back(br);
  border.push_back(bl);
  Polygon::Polygon_2 hole1;
  hole1.push_back(Point((xt+xmax)/2, (ymin+ymid)/2));
  hole1.push_back(Point(((xt+xmax)/2), ymid));
  hole1.push_back(Point(xt+(xmax-xt)/4, (ymin+ymid)/2));
  hole1.push_back(Point((xt+xmax)/2, (ymin+ymid)/2));
  Polygon::Polygon_2 hole2;
  hole2.push_back(Point((xt+xmin)/2, (ymin+ymid)/2));
  hole2.push_back(Point(((xt+xmin)/2), ymid));
  hole2.push_back(Point(xmin+(xt-xmin)/4, (ymin+ymid)/2));
  hole2.push_back(Point((xt+xmin)/2, (ymin+ymid)/2));
  Polygon::Holes_container holes;
  holes.push_back(hole1);
  holes.push_back(hole2);
  return Polygon(border, holes.begin(), holes.end());
}
MultiPoint generate_multipoint()
{
  MultiPoint mp;
  mp.push_back(generate_point(0,15,0,15));
  mp.push_back(generate_point(0,15,0,15));
  return mp;
}
MultiLinestring generate_multilinestring()
{
  Linestring l1=generate_linestring(),
      l2=generate_linestring();
  MultiLinestring mls;
  mls.push_back(l1);
  mls.push_back(l2);
  return mls;
}
MultiPolygon generate_multipolygon()
{
  Polygon p1=generate_polygon(),
      p2=generate_polygon();
  MultiPolygon polies;
  polies.push_back(p1);
  polies.push_back(p2);
  return polies;
}

int main()
{
  
  srand( unsigned(time(NULL) ));
  Point p = generate_point(0,6,0,6);
  {
    std::ofstream os("test.wkt");
    os.precision(17);
    CGAL::write_WKT(os, p);
    os.close();
  }
  Point test_p;
  {
    std::ifstream is("test.wkt");
    CGAL::read_point_WKT(is, test_p);
    is.close();
  }
  CGAL_assertion(p == test_p);
  
  Linestring ls = generate_linestring();
  {
    std::ofstream os("test.wkt");
    os.precision(17);
    CGAL::write_WKT(os, ls);
    os.close();
  }
  Linestring test_ls;
  {
    std::ifstream is("test.wkt");
    CGAL::read_linestring_WKT(is, test_ls);
    is.close();
  }
  CGAL_assertion(ls.range == test_ls.range);
  
  
  Polygon poly = generate_polygon();
  {
    std::ofstream os("test.wkt");
    os.precision(17);
    CGAL::write_WKT(os, poly);
    os.close();
  }
  Polygon test_poly;
  {
    std::ifstream is("test.wkt");
    CGAL::read_polygon_WKT(is, test_poly);
    is.close();
  }
  CGAL_assertion(poly == test_poly);
  
  MultiPoint pees = generate_multipoint();
  {
    std::ofstream os("test.wkt");
    os.precision(17);
    CGAL::write_WKT(os, pees);
    os.close();
  }
  MultiPoint test_pees;
  {
    std::ifstream is("test.wkt");
    CGAL::read_multipoint_WKT(is, test_pees);
    is.close();
  }
  CGAL_assertion(pees.range == test_pees.range);
  
  MultiLinestring mls = generate_multilinestring();
  {
    std::ofstream os("test.wkt");
    os.precision(17);
    CGAL::write_WKT(os, mls);
    os.close();
  }
  MultiLinestring test_mls;
  {
    std::ifstream is("test.wkt");
    CGAL::read_multilinestring_WKT(is, test_mls);
    is.close();
  }
  bool ok = true;
  for(size_t i=0; i<mls.size(); ++i)
    ok &= mls[i].range == test_mls[i].range;
  CGAL_assertion(ok);
  
  MultiPolygon polies = generate_multipolygon();
  {
    std::ofstream os("test.wkt");
    os.precision(17);
    CGAL::write_WKT(os, polies);
    os.close();
  }
  MultiPolygon test_polies;
  {
    std::ifstream is("test.wkt");
    CGAL::read_multipolygon_WKT(is, test_polies);
    is.close();
  }
  CGAL_assertion(polies.range == test_polies.range);
  
  std::cout<<"WKT writing test passed."<<std::endl;
  return 0;
}
