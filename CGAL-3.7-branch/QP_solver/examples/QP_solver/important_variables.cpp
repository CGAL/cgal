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
  // convex hull: 4-gon spanned by {(1,0), (4,1), (4,4), (2,3)}
  points.push_back (Point_d (1, 0)); // point 0
  points.push_back (Point_d (4, 1)); // point 1
  points.push_back (Point_d (4, 4)); // point 2
  points.push_back (Point_d (2, 3)); // point 3
 
  // test all 25 integer points in [0,4]^2
  for (int i=0; i<=4; ++i)
    for (int j=0; j<=4; ++j) {
      Point_d p (i, j);
      Solution s = solve_convex_hull_containment_lp
	(p, points.begin(), points.end(), CGAL::MP_Float());
      std::cout << p;
      if (s.is_infeasible())
	std::cout << " is not in the convex hull\n";
      else {
	assert (s.is_optimal());
	std::cout << " is a convex combination of the points ";
	Solution::Index_iterator it = s.basic_variable_indices_begin();
	Solution::Index_iterator end = s.basic_variable_indices_end();
	for (; it != end; ++it) std::cout << *it << " ";
	std::cout << std::endl;
      }
    }
  return 0;
}
