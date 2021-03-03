#include <iostream>

#include <CGAL/Epick_d.h>
#include <CGAL/Orthtree.h>
#include <CGAL/Orthtree_traits_d.h>
#include <CGAL/Random.h>

// Type Declarations
typedef CGAL::Dimension_tag<4> Dimension;
typedef CGAL::Epick_d<Dimension> Kernel;
typedef Kernel::Point_d Point_d;
typedef std::vector<Point_d> Point_vector;

typedef CGAL::Orthtree_traits_d<Kernel, Dimension> Traits;
typedef CGAL::Orthtree<Traits, Point_vector> Orthtree;

int main()
{
  CGAL::Random r;

  Point_vector points_dd;
  for (std::size_t i = 0; i < 5; ++ i)
  {
    std::array<double, Dimension::value> init;
    for (double& v : init)
      v = r.get_double(-1., 1.);
    points_dd.emplace_back (init.begin(), init.end());
  }

  Orthtree orthtree(points_dd);
  orthtree.refine(10, 5);

  return EXIT_SUCCESS;
}
