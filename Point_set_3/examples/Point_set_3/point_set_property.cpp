#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Point_set_3.h>
#include <CGAL/Random.h>

#include <fstream>
#include <limits>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::FT FT;
typedef Kernel::Point_3 Point;
typedef Kernel::Vector_3 Vector;
typedef std::array<unsigned char, 3> Color;

typedef CGAL::Point_set_3<Point> Point_set;
typedef Point_set::Property_map<Color> Color_map;
typedef Point_set::Property_map<FT> FT_map;

void print_point_set (const Point_set& point_set)
{
  Color_map color;
  boost::tie (color, boost::tuples::ignore) = point_set.property_map<Color>("color");
  FT_map intensity;
  boost::tie (intensity, boost::tuples::ignore) =  point_set.property_map<FT>("intensity");

  std::cerr << "Content of point set:" << std::endl;
  for (Point_set::const_iterator it = point_set.begin();
       it != point_set.end(); ++ it)
    std::cerr << "* Point " << point_set.point(*it) // or point_set[it]
              << " with color [" << (int)(color[*it][0])
              << " " << (int)(color[*it][1])
              << " " << (int)(color[*it][2])
              << "] and intensity " << intensity[*it]
              << std::endl;
}


int main (int, char**)
{
  Point_set point_set;

  Color black = {{ 0, 0, 0 }};
  bool success = false;
  Color_map color;

  boost::tie (color, success) = point_set.add_property_map<Color> ("color", black);
  assert (success);

  FT_map intensity;
  boost::tie (intensity, success) = point_set.add_property_map<FT> ("intensity", 0.);
  assert (success);

  point_set.reserve (10); // For memory optimization
  for (std::size_t i = 0; i < 10; ++ i)
    {
      Point_set::iterator it = point_set.insert (Point (double(i), double(i), double(i)));
      Color c = {{ (unsigned char)(CGAL::get_default_random().get_int(0, 255)),
                   (unsigned char)(CGAL::get_default_random().get_int(0, 255)),
                   (unsigned char)(CGAL::get_default_random().get_int(0, 255)) }};
      color[*it] = c;
      intensity[*it] = rand() / (double)(RAND_MAX);
    }

  print_point_set (point_set);

  // Remove points with intensity less than 0.5
  Point_set::iterator it = point_set.begin();
  while (it != point_set.end())
    {
      if (intensity[*it] < 0.5)
        point_set.remove(it);
      else
        ++ it;
    }

  print_point_set (point_set); // point set is filtered

  return 0;
}
