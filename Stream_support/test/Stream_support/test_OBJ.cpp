#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#include <CGAL/IO/OBJ.h>
#include <CGAL/IO/polygon_soup_io.h>
#include <CGAL/IO/OBJ.h>

#include <array>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

typedef CGAL::Simple_cartesian<double>                Kernel;
typedef Kernel::Point_3                               Point;
typedef std::vector<std::size_t>                      Face;

template <typename PointType, typename PolygonType>
void read(const std::filesystem::path& fname,
          std::size_t v, std::size_t f)
{
  std::cout << "Reading "<< fname << std::endl;
  std::ifstream input(fname);
  assert(input);

  std::cout << "Types: " << std::endl;
  std::cout << typeid(PointType).name() << std::endl;
  std::cout << typeid(PolygonType).name() << std::endl;

  std::cout << "Expecting " << v << " vertices and " << f << " triangles" << std::endl;

  std::vector<PointType> points;
  std::vector<PolygonType> faces;
  bool ok = CGAL::IO::read_OBJ(input, points, faces);

  std::cout << "Read " << points.size() << " vertices and " << faces.size() << " triangles" << std::endl;

  assert(ok);
  assert(points.size() == v);
  assert(faces.size() == f);
}

void test_types()
{
  // bunch of types to test
  typedef std::array<double, 3>                                       Point_type_1;
  typedef CGAL::Exact_predicates_exact_constructions_kernel::Point_3  Point_type_2;
  typedef std::basic_string<double>                                   Point_type_3;

  typedef std::array<int, 3>                                          Polygon_type_1;
  typedef std::vector<int>                                            Polygon_type_2;
  typedef std::basic_string<int>                                      Polygon_type_3;

  read<Point_type_1, Polygon_type_1>("data/90089.obj", 434, 864);
  read<Point_type_1, Polygon_type_2>("data/90089.obj", 434, 864);
  read<Point_type_1, Polygon_type_3>("data/90089.obj", 434, 864);
  read<Point_type_2, Polygon_type_1>(CGAL::data_file_path("meshes/pig.obj"), 468, 891);
  read<Point_type_2, Polygon_type_2>(CGAL::data_file_path("meshes/pig.obj"), 468, 891);
  read<Point_type_2, Polygon_type_3>(CGAL::data_file_path("meshes/pig.obj"), 468, 891);
  read<Point_type_3, Polygon_type_1>(CGAL::data_file_path("meshes/pig.obj"), 468, 891);
  read<Point_type_3, Polygon_type_2>(CGAL::data_file_path("meshes/pig.obj"), 468, 891);
  read<Point_type_3, Polygon_type_3>(CGAL::data_file_path("meshes/pig.obj"), 468, 891);
}

void test_out_of_range_indices()
{
  std::vector<Point> points;
  std::vector<Face> polygons;

  // valid: the last vertex, referenced by positive and by negative index
  std::istringstream valid("v 0 0 0\nv 1 0 0\nv 0 1 0\nf 1 2 3\nf -1 -2 -3\n");
  bool ok = CGAL::IO::read_OBJ(valid, points, polygons);
  assert(ok);
  assert(points.size() == 3 && polygons.size() == 2);

  // one past the last vertex, in a face
  points.clear();
  polygons.clear();
  std::istringstream face_past_end("v 0 0 0\nv 1 0 0\nv 0 1 0\nf 1 2 4\n");
  ok = CGAL::IO::read_OBJ(face_past_end, points, polygons);
  assert(!ok);

  // one past the last vertex, in a polyline
  points.clear();
  polygons.clear();
  std::vector<std::vector<std::size_t> > polylines;
  std::istringstream polyline_past_end("v 0 0 0\nv 1 0 0\nv 0 1 0\nl 1 4\n");
  ok = CGAL::IO::internal::read_OBJ(polyline_past_end, points, polygons, polylines);
  assert(!ok);

  // one before the first vertex, by negative index
  points.clear();
  polygons.clear();
  std::istringstream face_before_begin("v 0 0 0\nv 1 0 0\nv 0 1 0\nf -1 -2 -4\n");
  ok = CGAL::IO::read_OBJ(face_before_begin, points, polygons);
  assert(!ok);

  // negative indices count back from the vertices read so far, not from the end of the file
  points.clear();
  polygons.clear();
  std::istringstream relative_valid("v 0 0 0\nv 1 0 0\nv 0 1 0\nf -1 -2 -3\nv 1 1 1\n");
  ok = CGAL::IO::read_OBJ(relative_valid, points, polygons);
  assert(ok);
  assert(polygons.size() == 1 && polygons[0] == Face({2, 1, 0}));

  points.clear();
  polygons.clear();
  std::istringstream face_before_vertices("v 0 0 0\nf -1 -2 -3\nv 1 0 0\nv 0 1 0\n");
  ok = CGAL::IO::read_OBJ(face_before_vertices, points, polygons);
  assert(!ok);

  points.clear();
  polygons.clear();
  polylines.clear();
  std::istringstream polyline_before_vertices("v 0 0 0\nl -1 -2\nv 1 0 0\n");
  ok = CGAL::IO::internal::read_OBJ(polyline_before_vertices, points, polygons, polylines);
  assert(!ok);
}

int main(int argc, char** argv)
{
  test_out_of_range_indices();

  const char* obj_file = (argc > 1) ? argv[1] : "data/cube_quad.obj";

  std::vector<Point> points;
  std::vector<Face> polygons;

  bool ok = CGAL::IO::read_OBJ(obj_file, points, polygons, CGAL::parameters::verbose(true));
  assert(ok);
  std::cout << points.size() << " points and " << polygons.size() << " polygons" << std::endl;

  if(argc == 1)
    assert(points.size() == 8 && polygons.size() == 6);

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

  test_types();

  std::cout << "Done!" << std::endl;
  return EXIT_SUCCESS;
}
