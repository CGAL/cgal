#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#include <CGAL/array.h>
#include <CGAL/IO/STL_reader.h>

#include <cassert>
#include <fstream>
#include <iostream>
#include <vector>

template <typename Point_type, typename Polygon_type>
void read(const char* fname, std::size_t v, std::size_t f)
{
  std::cout << "Reading "<< fname << std::endl;
  std::ifstream input(fname, std::ios::in | std::ios::binary);

  std::cout << "Types: " << std::endl;
  std::cout << typeid(Point_type).name() << std::endl;
  std::cout << typeid(Polygon_type).name() << std::endl;

  std::cout << "Expecting " << v << " vertices and " << f << " triangles" << std::endl;

  std::vector<Point_type> points;
  std::vector<Polygon_type> faces;
  assert(CGAL::read_STL(input, points, faces, true));

  std::cout << "OFF version of file " << fname << std::endl;

  std::cout << "OFF\n" << points.size() << " " << faces.size()  << " 0" << std::endl;
  for(std::size_t i=0; i < points.size(); i++)
    std::cout << points[i][0] << " " << points[i][1] << " " << points[i][2]<< std::endl;

  for(std::size_t i=0; i < faces.size(); i++)
    std::cout << "3 " << faces[i][0] << " " << faces[i][1] << " " << faces[i][2] << std::endl;

  assert(points.size() == v);
  assert(faces.size() == f);
}

int main(int, char**)
{
  std::cout.precision(17);

  // bunch of types to test
  typedef CGAL::cpp11::array<double, 3>                               Point_type_1;
  typedef CGAL::Exact_predicates_exact_constructions_kernel::Point_3  Point_type_2;
  typedef std::basic_string<double>                                   Point_type_3;

  typedef CGAL::cpp11::array<int, 3>                                  Polygon_type_1;
  typedef std::vector<int>                                            Polygon_type_2;
  typedef std::basic_string<int>                                      Polygon_type_3;

  read<Point_type_1, Polygon_type_1>("data/cube.stl", 8, 12);
  read<Point_type_1, Polygon_type_2>("data/triangle.stl", 3, 1);

  read<Point_type_1, Polygon_type_3>("data/ascii-tetrahedron.stl", 4, 4);
  read<Point_type_2, Polygon_type_1>("data/binary-tetrahedron-nice-header.stl", 4, 4);
  read<Point_type_2, Polygon_type_2>("data/binary-tetrahedron-non-standard-header-1.stl", 4, 4);
  read<Point_type_2, Polygon_type_3>("data/binary-tetrahedron-non-standard-header-2.stl", 4, 4);
  read<Point_type_3, Polygon_type_1>("data/binary-tetrahedron-non-standard-header-3.stl", 4, 4);
  read<Point_type_3, Polygon_type_2>("data/binary-tetrahedron-non-standard-header-4.stl", 4, 4);
  read<Point_type_3, Polygon_type_3>("data/binary-tetrahedron-non-standard-header-5.stl", 4, 4);

  return EXIT_SUCCESS;
}
