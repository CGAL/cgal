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


  if (point_set.has_normals())
    std::cerr << "Point set has normals" << std::endl;
  else
    std::cerr << "Point set doesn't have normals" << std::endl;

  std::vector<Point_set::Item> indices;
  std::ifstream f (argc > 1 ? argv[1] : "data/data.pwn");

  
  Point_set point_set;
  point_set.add_normal_property();
  CGAL::read_xyz_points_and_normals(f,
                                    point_set.index_back_inserter(),
                                    point_set.point_push_pmap(),
                                    point_set.normal_push_pmap(),
                                    Kernel());
  
  CGAL::grid_simplify_point_set (point_set.begin (), point_set.end (),
                                 point_set.point_pmap(),
                                 0.1);

  typedef CGAL::cpp11::array<unsigned char, 3> Color;
  point_set.add_property<Color> ("color");

  for (std::size_t i = 0; i < (std::min)(std::size_t(5), point_set.size()); ++ i)
    std::cerr << "Item " << i << " = " << std::endl
              << "  * Point = " << point_set.point(i) << std::endl
              << "  * Point = " << point_set.normal(i) << std::endl
              << "  * Color = "
              << point_set.property<Color>("color", i)[0]
              << point_set.property<Color>("color", i)[1]
              << point_set.property<Color>("color", i)[2] << std::endl;

  
  f.close ();

  for (std::size_t i = 0; i < (std::min)(std::size_t(5), point_set.size()); ++ i)
    std::cerr << "Item " << i << " = " << std::endl
              << "  * Point = " << point_set[i] << std::endl
              << "  * Normal = " << point_set.normal(i) << std::endl;

  
  
  Point_set::iterator
    first_selected = CGAL::grid_simplify_point_set (point_set.begin (), point_set.end (),
                                                    &(point_set[0]),
                                                    0.1);

  std::cerr << std::endl;
  for (std::size_t i = 0; i < (std::min)(std::size_t(5), point_set.size()); ++ i)
    std::cerr << "Item " << i << " = " << std::endl
              << "  * Point = " << point_set[i] << std::endl
              << "  * Normal = " << point_set.normal(i) << std::endl;

  std::cerr << std::endl;
  std::size_t i = 0;
  for (Point_set::iterator it = first_selected;
       it != point_set.end() && i < 5; ++ i, ++ it)
    std::cerr << "Index " << *it << " = " << std::endl
              << "  * Point = " << point_set[*it] << std::endl
              << "  * Normal = " << point_set.normal(*it) << std::endl;

  if (!(point_set.are_indices_up_to_date()))
    std::cerr << "Indices not up to date" << std::endl;
  std::cerr << "Size = " << point_set.size() << std::endl;
  std::cerr << "Nb removed = " << std::distance (first_selected, point_set.end()) << std::endl;
  point_set.erase (first_selected, point_set.end());
  std::cerr << "Size = " << point_set.size() << std::endl;
  if (!(point_set.are_indices_up_to_date()))
    std::cerr << "Indices not up to date" << std::endl;

  for (std::size_t i = 0; i < (std::min)(std::size_t(5), point_set.size()); ++ i)
    std::cerr << "Item " << i << " = " << std::endl
              << "  * Point = " << point_set[i] << std::endl
              << "  * Normal = " << point_set.normal(i) << std::endl;

  std::ofstream fout ("out.xyz");
  CGAL::write_xyz_points (fout, point_set.begin(), point_set.end(), &(point_set[0]));
  fout.close();

  if (point_set.has_property<Color> ("color"))
    std::cerr << "Point set has colors" << std::endl;
  else
    std::cerr << "Point set doesn't have colors" << std::endl;
  
  point_set.add_property<Color> ("color");

  if (point_set.has_property<Color> ("color"))
    std::cerr << "Point set has colors" << std::endl;
  else
    std::cerr << "Point set doesn't have colors" << std::endl;

  for (std::size_t i = 0; i < (std::min)(std::size_t(5), point_set.size()); ++ i)
    std::cerr << "Item " << i << " = " << std::endl
              << "  * Point = " << point_set[i] << std::endl
              << "  * Color = " << point_set.property<Color>("color", i)[0] << std::endl;
  
  
  return 0;
}
