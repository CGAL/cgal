#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#include <CGAL/IO/GOCAD.h>
#include <CGAL/IO/polygon_soup_io.h>

#include <iostream>
#include <fstream>
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
  bool ok = CGAL::IO::read_GOCAD(input, points, faces);

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

  read<Point_type_1, Polygon_type_1>("data/2016206_MHT_surface.ts", 12491, 24191);
  read<Point_type_1, Polygon_type_2>("data/2016206_MHT_surface.ts", 12491, 24191);
  read<Point_type_1, Polygon_type_3>("data/2016206_MHT_surface.ts", 12491, 24191);
  read<Point_type_2, Polygon_type_1>("data/2016206_MHT_surface.ts", 12491, 24191);
  read<Point_type_2, Polygon_type_2>("data/2016206_MHT_surface.ts", 12491, 24191);
  read<Point_type_2, Polygon_type_3>("data/2016206_MHT_surface.ts", 12491, 24191);
  read<Point_type_3, Polygon_type_1>("data/2016206_MHT_surface.ts", 12491, 24191);
  read<Point_type_3, Polygon_type_2>("data/2016206_MHT_surface.ts", 12491, 24191);
  read<Point_type_3, Polygon_type_3>("data/2016206_MHT_surface.ts", 12491, 24191);
}

int main(int argc, char** argv)
{
  const char* gocad_file = (argc > 1) ? argv[1] : "data/2016206_MHT_surface.ts";

  std::vector<Point> points;
  std::vector<Face> polygons;

  bool ok = CGAL::IO::read_GOCAD(gocad_file, points, polygons);
  assert(ok);
  std::cout << points.size() << " points and " << polygons.size() << " polygons" << std::endl;

  if(argc == 0)
    assert(points.size() == 12491 && polygons.size() == 24191);

  points.clear();
  polygons.clear();
  std::string gocad_string(gocad_file);
  ok = CGAL::IO::read_GOCAD(gocad_string, points, polygons);
  assert(ok);

  points.clear();
  polygons.clear();
  std::ifstream is(gocad_file);
  ok = CGAL::IO::read_GOCAD(is, points, polygons);
  assert(ok);
  is.close();

  ok = CGAL::IO::write_GOCAD(gocad_file, points, polygons);
  assert(ok);

  std::ofstream os("tmp.ts");
  ok = CGAL::IO::write_GOCAD(os, points, polygons);
  assert(ok);
  os.close();

  ok = CGAL::IO::write_GOCAD("tmp.ts", points, polygons);
  assert(ok);

  const std::size_t ptn = points.size();
  const std::size_t pln = polygons.size();

  ok = CGAL::IO::write_polygon_soup("tmp.ts", points, polygons);
  assert(ok);

  points.clear();
  polygons.clear();

  ok = CGAL::IO::read_polygon_soup("tmp.ts", points, polygons);
  assert(ok);

  assert(points.size() == ptn);
  assert(polygons.size() == pln);

  test_types();

  std::cout << "Done!" << std::endl;
  return EXIT_SUCCESS;
}
