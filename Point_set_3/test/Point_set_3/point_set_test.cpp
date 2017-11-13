#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Point_set_3.h>
#include <CGAL/Point_set_3/Point_set_processing_3.h>
#include <CGAL/Point_set_3/IO.h>
#include <CGAL/grid_simplify_point_set.h>

#include <fstream>
#include <limits>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::FT FT;
typedef Kernel::Point_3 Point;
typedef Kernel::Vector_3 Vector;

typedef CGAL::Point_set_3<Point> Point_set;
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

  test (!(point_set.has_normal_map()), "point set shouldn't have normals.");
  point_set.add_normal_map();
  test (point_set.has_normal_map(), "point set should have normals.");

  std::ifstream f ("data/oni.pwn");
  CGAL::read_xyz_point_set(f, point_set);
  
  f.close ();

  Point_set::iterator
    first_to_remove = CGAL::grid_simplify_point_set (point_set.begin(),
                                                     point_set.end(),
                                                     point_set.point_map(),
                                                     0.1);

  std::size_t size = point_set.size ();
  point_set.remove_from (first_to_remove);
  test ((point_set.size() + point_set.garbage_size() == size), "sizes before and after removal do not match.");
 
  Point_set::Point_range
    range = point_set.points();

  {
    Point_set::const_iterator psit = point_set.begin();
    bool range_okay = true;
    for (Point_set::Point_range::const_iterator it = range.begin(); it != range.end(); ++ it)
      {
        if (*it != point_set.point (*psit))
          {
            range_okay = false;
            break;
          }
        ++ psit;
      }
    test (range_okay, "range access does not follow property map based access.");
  }

  test (point_set.has_garbage(), "point set should have garbage.");
  point_set.collect_garbage();
  test (!(point_set.has_garbage()), "point set shouldn't have garbage.");

  point_set.remove (*(point_set.begin()));
  
  test (point_set.has_garbage(), "point set should have garbage.");
  point_set.collect_garbage();
  test (!(point_set.has_garbage()), "point set shouldn't have garbage.");
  
  test (!(point_set.has_property_map<Color> ("color")), "point set shouldn't have colors.");
  Point_set::Property_map<Color> color_prop;
  bool garbage;
  boost::tie (color_prop, garbage) = point_set.add_property_map ("color", Color());
  test (point_set.has_property_map<Color> ("color"), "point set should have colors.");

  for (Point_set::iterator it = point_set.begin(); it != point_set.end(); ++ it)
    {
      Color c = {{ static_cast<unsigned char>(rand() % 255),
                   static_cast<unsigned char>(rand() % 255),
                   static_cast<unsigned char>(rand() % 255) }};
      put (color_prop, *it, c);
      test ((get (color_prop, *it) == c), "recovered color is incorrect.");
    }

  Point_set::Property_map<Color> color_prop_2;
  boost::tie (color_prop_2, garbage) = point_set.property_map<Color>("color");
  test ((color_prop_2 == color_prop), "color property not recovered correctly.");
  
  point_set.remove_normal_map ();
  test (!(point_set.has_normal_map()), "point set shouldn't have normals.");
  
  test (point_set.has_property_map<Color> ("color"), "point set should have colors.");
  point_set.remove_property_map<Color> (color_prop);
  test (!(point_set.has_property_map<Color> ("color")), "point set shouldn't have colors.");

  point_set.add_property_map<int> ("label", 0);
  point_set.add_property_map<double> ("intensity", 0.0);

  test (point_set.base().n_properties() == 4, "point set should have 4 properties.");
  
  Point p_before = *(point_set.points().begin());
  point_set.clear_properties();

  test (point_set.base().n_properties() == 2, "point set should have 2 properties.");
  test (!(point_set.has_property_map<int>("label")), "point set shouldn' have labels.");
  test (!(point_set.has_property_map<double>("intensity")), "point set shouldn' have intensity.");
  test (!(point_set.empty()), "point set shouldn' be empty.");

  Point p_after = *(point_set.points().begin());

  test (p_before == p_after, "points should not change when clearing properties.");

  std::cerr << nb_success << "/" << nb_test << " test(s) succeeded." << std::endl;
  
  return EXIT_SUCCESS;
}
