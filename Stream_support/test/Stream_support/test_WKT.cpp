#include <iostream>

#if BOOST_VERSION >= 105600 && (! defined(BOOST_GCC) || BOOST_GCC >= 40500)
#include <CGAL/IO/WKT.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <boost/config.hpp>
#include <boost/version.hpp>

#include <fstream>
#include <vector>

typedef CGAL::Exact_predicates_inexact_constructions_kernel           Kernel;

typedef CGAL::Point_2<Kernel>                                         Point;
typedef std::vector<Point>                                            Linestring;
typedef CGAL::Polygon_with_holes_2<Kernel>                            Poly;
typedef std::vector<Point>                                            MultiPoint;
typedef std::vector<Linestring>                                       MultiLinestring;
typedef std::vector<Poly>                                             MultiPolygon;

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
/// Read

bool test_read_WKT()
{
  Point p;
  {
    std::ifstream in("data/point.wkt");
    if(!CGAL::IO::read_point_WKT(in, p))
      return false;
    CGAL_assertion(p == Point(2,3));
  }
  {
    std::ifstream in("data/linestring.wkt");
    Linestring ls;
    if(!CGAL::IO::read_linestring_WKT(in, ls))
      return false;
    CGAL_assertion(ls.size() == 3);
  }
  {
    Poly poly;
    std::ifstream in("data/polygon.wkt");
    if(!CGAL::IO::read_polygon_WKT(in, poly))
      return false;
    CGAL_assertion(poly.outer_boundary().size() == 3);
  }
  {
    MultiPoint pees;
    std::ifstream in("data/multipoint.wkt");
    if(!CGAL::IO::read_multi_point_WKT(in, pees))
      return false;
    CGAL_assertion(pees.size() == 4);
  }
  {
    std::ifstream in("data/multilinestring.wkt");
    MultiLinestring mls;
    if(!CGAL::IO::read_multi_linestring_WKT(in, mls))
      return false;
    CGAL_assertion(mls.size() == 2);
  }
  {
    MultiPolygon polies;
    std::ifstream in("data/multipolygon.wkt");
    if(!CGAL::IO::read_multi_polygon_WKT(in, polies))
      return false;
    CGAL_assertion(polies.size() == 2);
  }

  std::cout << "WKT reading test passed." << std::endl;
  return true;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
// Write

double fRand(double fMin, double fMax)
{
  double f = (double)rand() / RAND_MAX;
  return fMin + f * (fMax - fMin);
}

Point generate_point(double xmin, double xmax,
                     double ymin, double ymax)
{
  double x = fRand(xmin, xmax),
         y = fRand(ymin, ymax);
  return Point(x,y);
}

Linestring generate_linestring()
{
  Linestring ls;
  ls.push_back(generate_point(0,15,0,15));
  ls.push_back(generate_point(0,15,0,15));
  return ls;
}

Poly generate_polygon()
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

  Poly::Polygon_2 border;
  border.push_back(bl);
  border.push_back(t);
  border.push_back(br);

  Poly::Polygon_2 hole1;
  hole1.emplace_back((xt+xmax)/2, (ymin+ymid)/2);
  hole1.emplace_back((xt+xmax)/2, ymid);
  hole1.emplace_back(xt+(xmax-xt)/4, (ymin+ymid)/2);

  Poly::Polygon_2 hole2;
  hole2.emplace_back((xt+xmin)/2, (ymin+ymid)/2);
  hole2.emplace_back((xt+xmin)/2, ymid);
  hole2.emplace_back(xmin+(xt-xmin)/4, (ymin+ymid)/2);

  Poly::Holes_container holes;
  holes.push_back(hole1);
  holes.push_back(hole2);

  return Poly(border, holes.begin(), holes.end());
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
  Linestring l1 = generate_linestring(),
             l2 = generate_linestring();

  MultiLinestring mls;
  mls.push_back(l1);
  mls.push_back(l2);

  return mls;
}

MultiPolygon generate_multipolygon()
{
  Poly p1 = generate_polygon(),
       p2 = generate_polygon();

  MultiPolygon polies;
  polies.push_back(p1);
  polies.push_back(p2);

  return polies;
}

bool test_write_WKT()
{
  srand( unsigned(time(NULL) ));

  Point p = generate_point(0,6,0,6);
  {
    std::ofstream os("test.wkt");
    os.precision(17);
    CGAL::IO::write_point_WKT(os, p);
    os.close();
  }
  Point test_p;
  {
    std::ifstream is("test.wkt");
    CGAL::IO::read_point_WKT(is, test_p);
    is.close();
  }
  CGAL_assertion(p == test_p);

  Linestring ls = generate_linestring();
  {
    std::ofstream os("test.wkt");
    os.precision(17);
    CGAL::IO::write_linestring_WKT(os, ls);
    os.close();
  }
  Linestring test_ls;
  {
    std::ifstream is("test.wkt");
    CGAL::IO::read_linestring_WKT(is, test_ls);
    is.close();
  }
  CGAL_assertion(ls == test_ls);

  Poly poly = generate_polygon();
  {
    std::ofstream os("test.wkt");
    os.precision(17);
    CGAL::IO::write_polygon_WKT(os, poly);
    os.close();
  }
  Poly test_poly;
  {
    std::ifstream is("test.wkt");
    CGAL::IO::read_polygon_WKT(is, test_poly);
    is.close();
  }

  CGAL_assertion(poly == test_poly);

  MultiPoint pees = generate_multipoint();
  {
    std::ofstream os("test.wkt");
    os.precision(17);
    CGAL::IO::write_multi_point_WKT(os, pees);
    os.close();
  }
  MultiPoint test_pees;
  {
    std::ifstream is("test.wkt");
    CGAL::IO::read_multi_point_WKT(is, test_pees);
    is.close();
  }
  CGAL_assertion(pees== test_pees);

  MultiLinestring mls = generate_multilinestring();
  {
    std::ofstream os("test.wkt");
    os.precision(17);
    CGAL::IO::write_multi_linestring_WKT(os, mls);
    os.close();
  }
  MultiLinestring test_mls;
  {
    std::ifstream is("test.wkt");
    CGAL::IO::read_multi_linestring_WKT(is, test_mls);
    is.close();
  }
  bool ok = true;
  for(size_t i=0; i<mls.size(); ++i)
    ok &= mls[i] == test_mls[i];
  CGAL_assertion(ok);

  MultiPolygon polies = generate_multipolygon();
  {
    std::ofstream os("test.wkt");
    os.precision(17);
    CGAL::IO::write_multi_polygon_WKT(os, polies);
    os.close();
  }
  MultiPolygon test_polies;
  {
    std::ifstream is("test.wkt");
    CGAL::IO::read_multi_polygon_WKT(is, test_polies);
    is.close();
  }
  CGAL_assertion(polies == test_polies);

  std::cout << "WKT writing test passed." << std::endl;
  return true;
}

int main()
{
  bool ok = test_read_WKT();
  assert(ok);
  ok = test_write_WKT();
  assert(ok);

  return EXIT_SUCCESS;
}
#else
int main(int, char**)
{
  return EXIT_SUCCESS;
}
#endif
