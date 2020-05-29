#include <CGAL/Simple_cartesian.h>

#include <CGAL/IO/GOCAD.h>
#include <CGAL/IO/polygon_soup_io.h>

#include <iostream>
#include <fstream>
#include <vector>

typedef CGAL::Simple_cartesian<double>                Kernel;
typedef Kernel::Point_3                               Point;
typedef std::vector<std::size_t>                      Face;

int main(int argc, char** argv)
{
  const char* ply_file = (argc > 1) ? argv[1] : "data/colored_tetra.ply";

  std::vector<Point> points;
  std::vector<Face> polygons;

  bool ok = CGAL::read_PLY(ply_file, points, polygons);
  assert(ok);
  std::cout << points.size() << " points and " << polygons.size() << " polygons" << std::endl;

  if(argc == 0)
    assert(points.size() == 4 && polygons.size() == 4);

  points.clear();
  polygons.clear();
  std::string ply_string(ply_file);
  ok = CGAL::read_PLY(ply_string, points, polygons);
  assert(ok);

  points.clear();
  polygons.clear();
  std::ifstream is(ply_file);
  ok = CGAL::read_PLY(is, points, polygons);
  assert(ok);
  is.close();

  ok = CGAL::write_PLY("tmp.ply", points, polygons);
  assert(ok);

  ok = CGAL::write_polygon_soup("tmp.ply", points, polygons);
  assert(ok);

  std::ofstream os("tmp.ply");
  CGAL::set_binary_mode(os);
  ok = CGAL::write_PLY(os, points, polygons);
  assert(ok);
  os.close();

  std::vector<Point> pts_backup = points;
  std::vector<Face> pls_backup = polygons;
  points.clear();
  polygons.clear();

  ok = CGAL::read_polygon_soup("tmp.ply", points, polygons);
  assert(ok);

  assert(points.size() == pts_backup.size());
  for(std::size_t i=0; i<points.size(); ++i)
    assert(CGAL::squared_distance(points[i], pts_backup[i]) < 1e-6);
  assert(polygons == pls_backup);

  std::cout << "Done!" << std::endl;
  return EXIT_SUCCESS;
}
