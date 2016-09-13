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


int main (int argc, char** argv)
{
  std::vector<Point_set::Item> indices;
  std::ifstream f (argc > 1 ? argv[1] : "data/data.pwn");

  Point_set point_set;
  point_set.add_normal_property();
  CGAL::read_xyz_points_and_normals(f,
                                    point_set.index_back_inserter(),
                                    point_set.point_push_pmap(),
                                    point_set.normal_push_pmap(),
                                    Kernel());
  std::cerr << point_set.size() << " point(s) read" << std::endl;
  
  Point_set::iterator first_to_remove
    = CGAL::grid_simplify_point_set (point_set.begin (), point_set.end (),
                                     point_set.point_pmap(),
                                     0.1);
  point_set.remove_from (first_to_remove);
  std::cerr << point_set.size() << " point(s) remaining after simplification" << std::endl;
  std::cerr << point_set.removed_size() << " removed point(s) remaining in memory" << std::endl;

  point_set.collect_garbage();
  
  std::cerr << point_set.removed_size() << " removed point(s) remaining after garbage collection" << std::endl;
  
  typename Point_set::Property_map<Color>::type color_prop;
  bool garbage;
  boost::tie (color_prop, garbage) = point_set.add_property ("color", Color());
  
  return 0;
}
