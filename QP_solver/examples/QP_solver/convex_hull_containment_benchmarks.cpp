// Example: assess the solver performance under any of the available 
// pricing strategies, in the convex-hull-containment problem
// NOTE: in order to see meaningful results, compile with -DNDEBUG
#include <vector>
#include <CGAL/Cartesian_d.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Random.h>
#include <CGAL/Timer.h>
#include "solve_convex_hull_containment_lp3.h"

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

int main()
{
  const int d = 10;       // change this in order to experiment
  const int n = 100000;   // change this in order to experiment

  // generate n random d-dimensional points in [0,1]^d
  CGAL::Random rd;
  std::vector<Point_d> points;
  for (int j =0; j<n; ++j) {
    std::vector<double> coords;
    for (int i=0; i<d; ++i) 
      coords.push_back(rd.get_double());
    points.push_back (Point_d (d, coords.begin(), coords.end()));
  }
  
  // benchmark all pricing strategies in turn
  CGAL::Quadratic_program_pricing_strategy strategy[] = {
    CGAL::QP_CHOOSE_DEFAULT,              // QP_PARTIAL_FILTERED_DANTZIG
    CGAL::QP_DANTZIG,                     // Dantzig's pivot rule...
    CGAL::QP_PARTIAL_DANTZIG,             // ... with partial pricing
    CGAL::QP_BLAND,                       // Bland's pivot rule
    CGAL::QP_FILTERED_DANTZIG,            // Dantzig's filtered pivot rule...
    CGAL::QP_PARTIAL_FILTERED_DANTZIG     // ... with partial pricing
  };
  
  CGAL::Timer t;
  for (int i=0; i<6; ++i) {
    // test strategy i
    CGAL::Quadratic_program_options options;
    options.set_pricing_strategy (strategy[i]);
    t.reset(); t.start();
    // is origin in convex hull of the points? (most likely, not)
    solve_convex_hull_containment_lp 
      (Point_d (d, CGAL::ORIGIN), points.begin(), points.end(), 
       ET(0), options);
    t.stop();
    std::cout << "Time (s) = " << t.time() << std::endl;
  }

  return 0;
}
