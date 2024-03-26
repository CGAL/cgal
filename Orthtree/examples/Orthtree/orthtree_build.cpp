#include <iostream>

#include <CGAL/Epick_d.h>
#include <CGAL/Cartesian_d.h>
#include <CGAL/Orthtree.h>
#include <CGAL/Orthtree_traits_point.h>
#include <CGAL/Random.h>

// Type Declarations
const int dimension = 4;
using Kernel = CGAL::Epick_d<CGAL::Dimension_tag<dimension> >;
using Point_d = Kernel::Point_d;
using Point_vector = std::vector<Point_d>;
using Traits = CGAL::Orthtree_traits_point<Kernel, Point_vector>;
using Orthtree = CGAL::Orthtree<Traits>;

int main()
{
  CGAL::Random r;

  Point_vector points_dd;
  for (std::size_t i = 0; i < 20; ++ i)
  {
    std::array<double, dimension> init{};
    for (double& v : init)
      v = r.get_double(-1., 1.);
    points_dd.emplace_back (init.begin(), init.end());
  }

  Orthtree orthtree(points_dd);
  orthtree.refine(10, 5);

  std::cout << orthtree << std::endl;

  return EXIT_SUCCESS;
}
