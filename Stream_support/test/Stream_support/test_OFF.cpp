#include <CGAL/Simple_cartesian.h>

#include <CGAL/IO/OFF.h>
#include <CGAL/IO/polygon_soup_io.h>

#include <fstream>
#include <iostream>
#include <vector>

typedef CGAL::Simple_cartesian<double>                Kernel;
typedef Kernel::Point_3                               Point;
typedef std::vector<std::size_t>                      Face;

int main(int argc, char** argv)
{
  const char* off_file = (argc > 1) ? argv[1] : "data/cube.off";

  std::vector<Point> points;
  std::vector<Face> polygons;

  bool ok = CGAL::IO::read_OFF(off_file, points, polygons);
  assert(ok);
  std::cout << points.size() << " points and " << polygons.size() << " polygons" << std::endl;

  if(argc == 0)
    assert(points.size() == 8 && polygons.size() == 12);

  points.clear();
  polygons.clear();
  std::string off_string(off_file);
  ok = CGAL::IO::read_OFF(off_string, points, polygons);
  assert(ok);

  points.clear();
  polygons.clear();
  std::ifstream is(off_file);
  ok = CGAL::IO::read_OFF(is, points, polygons);
  assert(ok);
  is.close();
  points.clear();
  polygons.clear();

  is.open(off_file, std::ios::binary);
  ok = CGAL::IO::read_OFF(is, points, polygons);
  assert(ok);
  is.close();

  std::ofstream os("tmp.off", std::ios::binary);
  ok = CGAL::IO::write_OFF(os, points, polygons);
  assert(ok);
  os.close();

  ok = CGAL::IO::write_OFF("tmp.off", points, polygons);
  assert(ok);

  std::vector<Point> pts_backup = points;
  std::vector<Face> pls_backup = polygons;

  ok = CGAL::IO::write_polygon_soup("tmp.off", points, polygons);
  assert(ok);

  points.clear();
  polygons.clear();

  ok = CGAL::IO::read_polygon_soup("tmp.off", points, polygons);
  assert(ok);

  assert(points.size() == pts_backup.size());
  for(std::size_t i=0; i<points.size(); ++i)
    assert(CGAL::squared_distance(points[i], pts_backup[i]) < 1e-6);
  assert(polygons == pls_backup);

  std::cout << "Done!" << std::endl;
  return EXIT_SUCCESS;
}
