// Example: find the points that contribute to a convex combination
#include <cassert>
#include <vector>
#include <CGAL/Cartesian_d.h>
#include <CGAL/MP_Float.h>
#include "solve_convex_hull_containment_lp2.h"

typedef CGAL::Cartesian_d<double> Kernel_d;
typedef Kernel_d::Point_d Point_d;
typedef CGAL::Quadratic_program_solution<CGAL::MP_Float> Solution;

int main()
{
  std::vector<Point_d> points;
  // convex hull: line spanned by {(0,0), (1,0),..., (9,0)}
  for (int j=0; j<10; ++j)
    points.push_back (Point_d (j, 0));

  for (double f=0.5; f<10; ++f) {
    Point_d p (f, 0.0);
    Solution s = solve_convex_hull_containment_lp
      (p, points.begin(), points.end(), CGAL::MP_Float());
    std::cout << p;
    if (s.status() == CGAL::QP_INFEASIBLE)
      std::cout << " is not in the convex hull." << std::endl;
    else {
      std::cout << " is a convex combination of the points ";
      Solution::Index_iterator it = s.basic_variable_indices_begin();
      Solution::Index_iterator end = s.basic_variable_indices_end();
      for (; it != end; ++it)
	std::cout << points[*it] << " ";
      std::cout << std::endl;
    }
  }

  return 0;
}
