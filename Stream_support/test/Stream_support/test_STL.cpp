#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#include <CGAL/IO/STL.h>
#include <CGAL/IO/polygon_soup_io.h>

#include <iostream>
#include <fstream>
#include <vector>

typedef CGAL::Exact_predicates_inexact_constructions_kernel   K;
typedef K::Point_3                                            Point;
typedef std::vector<std::size_t>                              Face;

template <typename Point_type, typename Polygon_type>
void read(const char* fname, std::size_t v, std::size_t f,
          bool is_binary = false, bool should_fail = false)
{
  std::cout << "Reading "<< fname << std::endl;
  std::ifstream input(fname, std::ios::in | std::ios::binary);

  std::cout << "Types: " << std::endl;
  std::cout << typeid(Point_type).name() << std::endl;
  std::cout << typeid(Polygon_type).name() << std::endl;

  std::cout << "Expecting " << v << " vertices and " << f << " triangles" << std::endl;

  std::vector<Point_type> points;
  std::vector<Polygon_type> faces;
  bool ok = CGAL::IO::read_STL(input, points, faces, CGAL::parameters::use_binary_mode(is_binary));
  assert(ok != should_fail);
  if(!should_fail)
  {
    std::cout << "OFF version of file " << fname << std::endl;

    std::cout << "OFF\n" << points.size() << " " << faces.size()  << " 0" << std::endl;
    for(std::size_t i=0; i < points.size(); i++)
      std::cout << points[i][0] << " " << points[i][1] << " " << points[i][2]<< std::endl;

    for(std::size_t i=0; i < faces.size(); i++)
      std::cout << "3 " << faces[i][0] << " " << faces[i][1] << " " << faces[i][2] << std::endl;

    assert(points.size() == v);
    assert(faces.size() == f);
  }
}

void further_tests()
{
  // bunch of types to test
  typedef std::array<double, 3>                                       Point_type_1;
  typedef CGAL::Exact_predicates_exact_constructions_kernel::Point_3  Point_type_2;
  typedef std::basic_string<double>                                   Point_type_3;

  typedef std::array<int, 3>                                          Polygon_type_1;
  typedef std::vector<int>                                            Polygon_type_2;
  typedef std::basic_string<int>                                      Polygon_type_3;

  read<Point_type_1, Polygon_type_1>("data/cube.stl", 8, 12, true);
  read<Point_type_1, Polygon_type_2>("data/triangle.stl", 3, 1);

  read<Point_type_1, Polygon_type_3>("data/ascii-tetrahedron.stl", 4, 4);
  read<Point_type_2, Polygon_type_1>("data/binary-tetrahedron-nice-header.stl", 4, 4, true);
  read<Point_type_2, Polygon_type_2>("data/binary-tetrahedron-non-standard-header-1.stl", 4, 4, true);
  read<Point_type_2, Polygon_type_3>("data/binary-tetrahedron-non-standard-header-2.stl", 4, 4, true);
  read<Point_type_3, Polygon_type_1>("data/binary-tetrahedron-non-standard-header-3.stl", 4, 4, true);
  read<Point_type_3, Polygon_type_2>("data/binary-tetrahedron-non-standard-header-4.stl", 4, 4, true);
  read<Point_type_3, Polygon_type_3>("data/binary-tetrahedron-non-standard-header-5.stl", 4, 4, true);
}

int main(int argc, char** argv)
{
  const char* stl_file = (argc > 1) ? argv[1] : "data/ascii-tetrahedron.stl";

  std::vector<Point> points;
  std::vector<Face> polygons;

  bool ok = CGAL::IO::read_STL(stl_file, points, polygons, CGAL::parameters::verbose(true));
  assert(ok);
  std::cout << points.size() << " points and " << polygons.size() << " polygons" << std::endl;

  if(argc == 0)
    assert(points.size() == 434 && polygons.size() == 864);

  points.clear();
  polygons.clear();
  std::string stl_string(stl_file);
  ok = CGAL::IO::read_STL(stl_string, points, polygons);
  assert(ok);

  points.clear();
  polygons.clear();
  std::ifstream is(stl_file);
  ok = CGAL::IO::read_STL(is, points, polygons);
  assert(ok);
  points.clear();
  polygons.clear();
  is.clear();
  is.seekg(0, is.beg);
  ok = CGAL::IO::read_STL(is, points, polygons, CGAL::parameters::use_binary_mode(false));
  assert(ok);
  is.close();

  ok = CGAL::IO::write_STL("tmp.stl", points, polygons);
  assert(ok);

  ok = CGAL::IO::write_polygon_soup("tmp.stl", points, polygons);
  assert(ok);

  std::ofstream os("tmp.stl");
  CGAL::IO::set_binary_mode(os);
  ok = CGAL::IO::write_STL(os, points, polygons);
  assert(ok);
  os.close();

  std::vector<Point> pts_backup = points;
  std::vector<Face> pls_backup = polygons;

  points.clear();
  polygons.clear();

  ok = CGAL::IO::read_polygon_soup("tmp.stl", points, polygons);
  assert(ok);

  assert(points.size() == pts_backup.size());
  for(std::size_t i=0; i<points.size(); ++i)
    assert(CGAL::squared_distance(points[i], pts_backup[i]) < 1e-6);
  assert(polygons == pls_backup);

  further_tests();

  std::cout << "Done!" << std::endl;
  return EXIT_SUCCESS;
}
