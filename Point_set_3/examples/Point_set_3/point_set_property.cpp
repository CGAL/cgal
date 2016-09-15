#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Point_set_3.h>

#include <fstream>
#include <limits>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::FT FT;
typedef Kernel::Point_3 Point;
typedef Kernel::Vector_3 Vector;
typedef CGAL::cpp11::array<unsigned char, 3> Color;

typedef CGAL::Point_set_3<Kernel> Point_set;
typedef Point_set::Property_map<Color>::type Color_prop;
typedef Point_set::Property_map<FT>::type Floating_prop;

void print_point_set (const Point_set& point_set)
{
  Color_prop color;
  boost::tie (color, boost::tuples::ignore) = point_set.get_property<Color>("color");
  Floating_prop intensity;
  boost::tie (intensity, boost::tuples::ignore) =  point_set.get_property<FT>("intensity");
  
  std::cerr << "Content of point set:" << std::endl;
  for (Point_set::const_iterator it = point_set.begin();
       it != point_set.end(); ++ it)
    std::cerr << "* Point " << point_set.point(it) // or point_set[it]
              << " with color [" << (int)(color[*it][0])
              << " " << (int)(color[*it][1])
              << " " << (int)(color[*it][2])
              << "] and intensity " << intensity[*it]
              << std::endl;
}


int main (int, char**)
{
  srand (time (NULL));
  
  Point_set point_set;

  Color black = {{ 0, 0, 0 }};
  bool okay = false;
  Color_prop color;

  boost::tie (color, okay) = point_set.add_property<Color> ("color", black);
  assert (okay);

  Floating_prop intensity;
  boost::tie (intensity, okay) = point_set.add_property<FT> ("intensity", 0.);
  assert (okay);

  point_set.reserve (10); // For memory optimization
  for (std::size_t i = 0; i < 10; ++ i)
    {
      Point_set::iterator it = point_set.push_back (Point (i, i, i));
      Color c = {{ (unsigned char)(rand() % 256),
                   (unsigned char)(rand() % 256),
                   (unsigned char)(rand() % 256) }};
      color[*it] = c;
      intensity[*it] = rand() / (double)(RAND_MAX);
    }

  print_point_set (point_set);

  std::random_shuffle (point_set.begin(), point_set.end());

  print_point_set (point_set); // point set is randomized

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
