#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Point_set_3.h>
#include <CGAL/IO/read_xyz_points.h>
#include <CGAL/IO/write_xyz_points.h>
#include <CGAL/grid_simplify_point_set.h>

#include <fstream>
#include <limits>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::FT FT;
typedef Kernel::Point_3 Point;
typedef Kernel::Vector_3 Vector;

typedef CGAL::Point_set_3<Kernel> Point_set;
typedef CGAL::cpp11::array<unsigned char, 3> Color;


void test (bool expr, const char* msg)
{
  static std::size_t nb_test = 0;
  ++ nb_test;
  if (!expr)
    std::cerr << "Error on test " << nb_test << ": " << msg << std::endl;
}


int main (int, char**)
{
  Point_set point_set;

  test (!(point_set.has_normals()), "point set shouldn't have normals.");
  point_set.add_normal_property();
  test (point_set.has_normals(), "point set should have normals.");

  std::vector<Point_set::Item> indices;
  std::ifstream f ("data/oni.pwn");
  CGAL::read_xyz_points_and_normals(f,
                                    point_set.index_back_inserter(),
                                    point_set.point_push_pmap(),
                                    point_set.normal_push_pmap(),
                                    Kernel());
  f.close ();

  Point_set::iterator
    first_to_remove = CGAL::grid_simplify_point_set (point_set.begin(),
                                                     point_set.end(),
                                                     &(point_set[0]),
                                                     0.1);

  std::size_t size = point_set.size ();
  point_set.remove_from (first_to_remove);
  test ((point_set.size() + point_set.removed_size() == size), "sizes before and after removal do not match.");
  
  test (point_set.has_garbage(), "point set should have garbage.");
  point_set.collect_garbage();
  test (!(point_set.has_garbage()), "point set shouldn't have garbage.");
  
  test (!(point_set.has_property<Color> ("color")), "point set shouldn't have colors.");
  point_set.add_property<Color> ("color");
  test (point_set.has_property<Color> ("color"), "point set should have colors.");


  for (std::size_t i = 0; i < point_set.size(); ++ i)
    {
      Color c = {{ static_cast<unsigned char>(rand() % 255),
                   static_cast<unsigned char>(rand() % 255),
                   static_cast<unsigned char>(rand() % 255) }};
      point_set.property<Color> ("color", i) = c;
      test ((point_set.property<Color> ("color", i) == c), "recovered color is incorrect.");
    }
  
  point_set.remove_normal_property ();
  test (!(point_set.has_normals()), "point set shouldn't have normals.");
  
  test (point_set.has_property<Color> ("color"), "point set should have colors.");
  point_set.remove_property<Color> ("color");
  test (!(point_set.has_property<Color> ("color")), "point set shouldn't have colors.");
 
  return 0;
}
