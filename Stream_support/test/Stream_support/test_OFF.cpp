#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#include <CGAL/IO/OFF.h>
#include <CGAL/IO/polygon_soup_io.h>

#include <array>
#include <fstream>
#include <iostream>
#include <vector>

typedef CGAL::Simple_cartesian<double>                Kernel;
typedef Kernel::Point_3                               Point;
typedef std::vector<std::size_t>                      Face;

template <typename PointType, typename PolygonType>
void read(const std::string& fname,
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
  bool ok = CGAL::IO::read_OFF(input, points, faces);

  std::cout << "Read " << points.size() << " vertices and " << faces.size() << " triangles" << std::endl;

  assert(ok);
  assert(points.size() == v);
  assert(faces.size() == f);
}

void test_types()
{
  // OFF currently uses Kernel_traits<Point_type>
  //
  // typedef std::array<double, 3>                                       Point_type_1;
  typedef CGAL::Exact_predicates_exact_constructions_kernel::Point_3  Point_type_2;
  // typedef std::basic_string<double>                                   Point_type_3;

  typedef std::array<int, 3>                                          Polygon_type_1;
  typedef std::vector<int>                                            Polygon_type_2;
  typedef std::basic_string<int>                                      Polygon_type_3;

  read<Point_type_2, Polygon_type_1>(CGAL::data_file_path("meshes/cactus.off"), 620, 1236);
  read<Point_type_2, Polygon_type_2>(CGAL::data_file_path("meshes/mpi_and_sphere.off"), 252, 372);
  read<Point_type_2, Polygon_type_3>(CGAL::data_file_path("meshes/blobby-shuffled.off"), 2027, 4050);
}

int main(int argc, char** argv)
{
  const std::string off_file = (argc > 1) ? argv[1] : CGAL::data_file_path("meshes/cube.off");

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

  test_types();

  std::cout << "Done!" << std::endl;
  return EXIT_SUCCESS;
}
