// Example: check whether a point is in the convex hull of other points
#include <cassert>
#include <vector>
#include <CGAL/Cartesian_d.h>
#include <CGAL/MP_Float.h>
#include "solve_convex_hull_containment_lp2.h"

// choose exact floating-point type
#ifdef CGAL_USE_GMP
#include <CGAL/Gmpzf.h>
typedef CGAL::Gmpzf ET;
#else
#include <CGAL/MP_Float.h>
typedef CGAL::MP_Float ET;
#endif

typedef CGAL::Cartesian_d<double> Kernel_d;
typedef Kernel_d::Point_d Point_d;

bool is_in_convex_hull (const Point_d& p,
			std::vector<Point_d>::const_iterator begin,
			std::vector<Point_d>::const_iterator end)
{
  CGAL::Quadratic_program_solution<ET> s =
    solve_convex_hull_containment_lp (p, begin, end,ET(0));
  return !s.is_infeasible();
}

int main()
{
  std::vector<Point_d> points;
  // convex hull: simplex spanned by {(0,0), (10,0), (0,10)}
  points.push_back (Point_d ( 0.0,  0.0));
  points.push_back (Point_d (10.0,  0.0));
  points.push_back (Point_d ( 0.0, 10.0));
  for (int i=0; i<=10; ++i)
    for (int j=0; j<=10; ++j) {
      // (i,j) is in the simplex iff i+j <= 10
      bool contained = is_in_convex_hull
	(Point_d (i, j), points.begin(), points.end());
      assert (contained == (i+j<=10));
    }

  return 0;
}
