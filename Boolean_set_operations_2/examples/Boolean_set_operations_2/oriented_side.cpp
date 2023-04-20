/*! \file oriented_side.cpp
 * Compute the oriented side of a point and a polygon with holes.
 */

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Boolean_set_operations_2.h>

typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
typedef Kernel::Point_2                                   Point_2;
typedef CGAL::Polygon_2<Kernel>                           Polygon_2;
typedef CGAL::Polygon_with_holes_2<Kernel>                Polygon_with_holes_2;

# define nice(os) ((os == CGAL::ON_ORIENTED_BOUNDARY) ? "on boundary" :  \
                   (os == CGAL::POSITIVE) ? "inside" : "outside")

int main() {
  Polygon_2 hole;
  hole.push_back(Point_2(1, 1));
  hole.push_back(Point_2(1, 2));
  hole.push_back(Point_2(2, 2));
  hole.push_back(Point_2(2, 1));

  Polygon_2 out;
  out.push_back(Point_2(0, 0));
  out.push_back(Point_2(3, 0));
  out.push_back(Point_2(3, 3));
  out.push_back(Point_2(0, 3));

  Polygon_with_holes_2 pwh(out, &hole, &hole+1);
  std::cout << pwh << std::endl;

  auto os = CGAL::oriented_side(Point_2(0, 0), pwh);
  std::cout << "(0,0) is : " << nice(os) << std::endl;
  os = CGAL::oriented_side(Point_2(0.5, 0.5), pwh);
  std::cout << "(0,0) is : " << nice(os) << std::endl;
  os = CGAL::oriented_side(Point_2(1, 1), pwh);
  std::cout << "(0,0) is : " << nice(os) << std::endl;
  os = CGAL::oriented_side(Point_2(2.5, 2.5), pwh);
  std::cout << "(0,0) is : " << nice(os) << std::endl;
  os = CGAL::oriented_side(Point_2(3, 3), pwh);
  std::cout << "(0,0) is : " << nice(os) << std::endl;
  os = CGAL::oriented_side(Point_2(3.5, 3.5), pwh);
  std::cout << "(0,0) is : " << nice(os) << std::endl;
  return 0;
}
