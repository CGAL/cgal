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

std::size_t nb_test = 0;
std::size_t nb_success = 0;

void test (bool expr, const char* msg)
{
  ++ nb_test;
  if (!expr)
    std::cerr << "Error on test " << nb_test << ": " << msg << std::endl;
  else
    ++ nb_success;
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
                                                     point_set.point_pmap(),
                                                     0.1);

  std::size_t size = point_set.size ();
  point_set.remove_from (first_to_remove);
  test ((point_set.size() + point_set.removed_size() == size), "sizes before and after removal do not match.");
  
  test (point_set.has_garbage(), "point set should have garbage.");
  point_set.collect_garbage();
  test (!(point_set.has_garbage()), "point set shouldn't have garbage.");
  
  test (!(point_set.has_property<Color> ("color")), "point set shouldn't have colors.");
  typename Point_set::Property_map<Color>::type color_prop;
  bool garbage;
  boost::tie (color_prop, garbage) = point_set.add_property ("color", Color());
  test (point_set.has_property<Color> ("color"), "point set should have colors.");

  for (Point_set::iterator it = point_set.begin(); it != point_set.end(); ++ it)
    {
      Color c = {{ static_cast<unsigned char>(rand() % 255),
                   static_cast<unsigned char>(rand() % 255),
                   static_cast<unsigned char>(rand() % 255) }};
      put (color_prop, *it, c);
      test ((get (color_prop, *it) == c), "recovered color is incorrect.");
    }

  typename Point_set::Property_map<Color>::type color_prop_2;
  boost::tie (color_prop_2, garbage) = point_set.property<Color>("color");
  test ((color_prop_2 == color_prop), "color property not recovered correctly.");
  
  point_set.remove_normal_property ();
  test (!(point_set.has_normals()), "point set shouldn't have normals.");
  
  test (point_set.has_property<Color> ("color"), "point set should have colors.");
  point_set.remove_property (color_prop);
  test (!(point_set.has_property<Color> ("color")), "point set shouldn't have colors.");

  std::cerr << nb_success << "/" << nb_test << " test(s) succeeded." << std::endl;
  
  return EXIT_SUCCESS;
}
