#include <iostream>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Quadtree.h>
#include <CGAL/Random.h>

// Type Declarations
using Kernel = CGAL::Simple_cartesian<double>;
using Point_2 = Kernel::Point_2;
using Point_vector = std::vector<Point_2>;
using Quadtree = CGAL::Quadtree<Kernel, Point_vector>;

int main()
{
  CGAL::Random r;
  Point_vector points_2d;
  for (std::size_t i = 0; i < 5; ++ i)
    points_2d.emplace_back(r.get_double(-1., 1.),
                           r.get_double(-1., 1.));

  Quadtree quadtree(points_2d);
  quadtree.refine(10, 5);

  return EXIT_SUCCESS;
}
