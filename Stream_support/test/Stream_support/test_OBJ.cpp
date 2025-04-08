#include <CGAL/Simple_cartesian.h>

#include <CGAL/IO/OBJ.h>
#include <CGAL/IO/polygon_soup_io.h>

#include <iostream>
#include <fstream>
#include <vector>

typedef CGAL::Simple_cartesian<double>                Kernel;
typedef Kernel::Point_3                               Point;
typedef std::vector<std::size_t>                      Face;
typedef Kernel::Point_2 Texture;
typedef Kernel::Vector_3 Normal;

int main(int argc, char** argv)
{
  const char* obj_file = (argc > 1) ? argv[1] : "data/90089.obj";

  std::vector<Point> points;
  std::vector<Face> polygons;
  std::vector<Normal> normals;
  std::vector<Texture> UVcoords;
  std::vector<int> normal_indices;
  std::vector<int> texture_indices;

  bool ok = CGAL::IO::read_OBJ(obj_file, points, polygons, CGAL::parameters::verbose(true));
  assert(ok);
  std::cout << points.size() << " points and " << polygons.size() << " polygons" << std::endl;

  if(argc == 1)
    assert(points.size() == 434 && polygons.size() == 864);

  points.clear();
  polygons.clear();

  // Check for internal read_OBJ with normal and texture coordinates reading capapbility
  std::ifstream input(obj_file);
  ok = CGAL::IO::internal::read_OBJ(input, points, polygons, std::back_inserter(normals), std::back_inserter(UVcoords),
                                    std::back_inserter(normal_indices), std::back_inserter(texture_indices), true);
  assert(ok);
  std::cout << points.size() << " points " << polygons.size() << " polygons " << normals.size() << " normals "
            << UVcoords.size() << " texture coords" << std::endl;

  if(argc == 1)
    assert(points.size() == 434 && polygons.size() == 864 && normals.size() == 242 && UVcoords.size() == 627);
  
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

  std::cout << "Done!" << std::endl;
  return EXIT_SUCCESS;
}
