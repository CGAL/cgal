#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Point_set_3.h>

#include <fstream>
#include <limits>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::FT FT;
typedef Kernel::Point_3 Point;
typedef Kernel::Vector_3 Vector;

typedef CGAL::Point_set_3<Point> Point_set;

void print_point_set (const Point_set& point_set)
{
  std::cerr << "Content of point set:" << std::endl;
  for (Point_set::const_iterator it = point_set.begin();
       it != point_set.end(); ++ it)
    std::cerr << "* Point " << point_set.point(*it) // or point_set[it]
              << " with normal " << point_set.normal(*it)
              << std::endl;
}

int main (int, char**)
{
  Point_set point_set;
  //  Point_set point_set2 ("visibility",Vector(0,0,0), "insensity", 1.0);
  //  std::cerr << point_set2.properties();
  
  // Add points
  point_set.insert (Point (0., 0., 0.));
  point_set.insert (Point (0., 0., 1.));
  point_set.insert (Point (0., 1., 0.));

  point_set.add_normal_map();

  print_point_set(point_set); // Normals have default values

  // Change normal values
  point_set.normal(0) = Vector (1., 0., 0.);
  point_set.normal(1) = Vector (0., 1., 0.);
  point_set.normal(2) = Vector (0., 0., 1.);

  // Add point + normal
  point_set.insert (Point (1., 2., 3.), Vector (4., 5., 6.));

  print_point_set(point_set); 

  // Add new item
  Point_set::iterator new_item = point_set.insert(Point (7., 8., 9.));
  point_set.normal(*new_item) = Vector (10., 11., 12.);

  print_point_set(point_set); // New item has default values

  point_set.remove (point_set.begin());

  print_point_set(point_set); // New item has default values
  
  // Display information
  std::cerr << point_set.properties();

  point_set.collect_garbage();
  
  // Display information (garbage should be gone)
  std::cerr << point_set.properties();

  return 0;
}
