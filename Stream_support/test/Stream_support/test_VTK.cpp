#include <CGAL/Simple_cartesian.h>

#include <CGAL/IO/VTK.h>

#include <iostream>
#include <fstream>
#include <vector>
#include <list>

typedef CGAL::Simple_cartesian<double>                Kernel;
typedef Kernel::Point_3                               Point;
typedef std::vector<long unsigned int>                Face;

int main()
{
  const char* vtk_file = "test_VTK.vtp";

  std::vector<Point> points =  { Point(0,0,0), Point(1,0,0), Point(0,1,0) };
  std::vector<Face> polygons = { Face{0,1,2} };

  bool ok = CGAL::IO::write_VTP(vtk_file, points, polygons, CGAL::parameters::verbose(true).use_binary_mode(false));
  assert(ok);

  points.clear();
  polygons.clear();
  ok = CGAL::IO::read_VTP(vtk_file, points, polygons, CGAL::parameters::verbose(true));
  assert(ok);
  std::cout << points.size() << " points and " << polygons.size() << " polygons" << std::endl;
#if 0
  if(argc == 1)
    assert(points.size() == 434 && polygons.size() == 864);

  points.clear();
  polygons.clear();
  std::string obj_string(obj_file);
  ok = CGAL::IO::read_OBJ(obj_string, points, polygons);
  assert(ok);

  points.clear();
  polygons.clear();

  std::ifstream is(obj_file);
  ok = CGAL::IO::read_OBJ(is, points, polygons);
  assert(ok);
  is.close();

  is.open(obj_file, std::ios::binary);
  ok = CGAL::IO::read_OBJ(is, points, polygons);
  assert(ok);
  is.close();

  std::ofstream os("tmp.obj");
  ok = CGAL::IO::write_OBJ(os, points, polygons);
  assert(ok);
  os.close();

  ok = CGAL::IO::write_OBJ("tmp.obj", points, polygons);
  assert(ok);

  std::vector<Point> pts_backup = points;
  std::vector<Face> pls_backup = polygons;

  ok = CGAL::IO::write_polygon_soup("tmp.obj", points, polygons);
  assert(ok);

  points.clear();
  polygons.clear();

  ok = CGAL::IO::read_polygon_soup("tmp.obj", points, polygons);
  assert(ok);

  assert(points.size() == pts_backup.size());
  for(std::size_t i=0; i<points.size(); ++i)
    assert(CGAL::squared_distance(points[i], pts_backup[i]) < 1e-6);
  assert(polygons == pls_backup);
#endif

  std::cout << "Done!" << std::endl;
  return EXIT_SUCCESS;
}
