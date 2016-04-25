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


int main (int, char**)
{
  Point_set point_set;

  assert (!(point_set.has_normals()));
  point_set.add_normal_property();
  assert (point_set.has_normals());

  std::vector<Point_set::Item> indices;
  std::ifstream f ("data/oni.pwn");
  CGAL::read_xyz_points_and_normals(f,
                                    point_set.index_back_inserter(),
                                    point_set.point_push_pmap(),
                                    point_set.normal_push_pmap(),
                                    Kernel());
  f.close ();

  Point_set::iterator
    first_selected = CGAL::grid_simplify_point_set (point_set.begin (), point_set.end (),
                                                    &(point_set[0]),
                                                    0.1);

  assert (!(point_set.are_indices_up_to_date()));
  point_set.erase (first_selected, point_set.end());
  assert (point_set.are_indices_up_to_date());

  assert (!(point_set.has_property<Color> ("color")));
  point_set.add_property<Color> ("color");
  assert (point_set.has_property<Color> ("color"));


  for (std::size_t i = 0; i < point_set.size(); ++ i)
    {
      Color c = {{ static_cast<unsigned char>(rand() % 255),
                   static_cast<unsigned char>(rand() % 255),
                   static_cast<unsigned char>(rand() % 255) }};
      point_set.property<Color> ("color", i) = c;
      assert (point_set.property<Color> ("color", i) == c);
    }
  
  point_set.remove_normal_property ();
  assert (!(point_set.has_normals()));
  
  assert (point_set.has_property<Color> ("color"));
  point_set.remove_property<Color> ("color");
  assert (!(point_set.has_property<Color> ("color")));
 
  return 0;
}
