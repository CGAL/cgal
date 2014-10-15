//#undef CGAL_LINKED_WITH_TBB // CJTODO TEMP

// Without TBB_USE_THREADING_TOOL Intel Inspector XE will report false positives in Intel TBB
// (http://software.intel.com/en-us/articles/compiler-settings-for-threading-error-analysis-in-intel-inspector-xe/)
#ifdef _DEBUG
# define TBB_USE_THREADING_TOOL
#endif

#include "test_utilities.h"

#include <CGAL/Epick_d.h>
#include <CGAL/Tangential_complex.h>
#include <CGAL/Random.h>
#include <CGAL/Kernel_traits.h>
#include <CGAL/Mesh_3/Profiling_tools.h>

#include <fstream>
#include <math.h>

#ifdef CGAL_LINKED_WITH_TBB
# include <tbb/task_scheduler_init.h>
#endif

#ifdef _DEBUG
  const int NUM_POINTS = 50;
#else
  const int NUM_POINTS = 500;
#endif

int main()
{
  const int INTRINSIC_DIMENSION = 4;
  const int AMBIENT_DIMENSION = 5;

  typedef CGAL::Epick_d<CGAL::Dimension_tag<AMBIENT_DIMENSION> > Kernel;
  typedef Kernel::Point_d                                        Point;
 
  //CGAL::default_random = CGAL::Random(0); // NO RANDOM

#ifdef CGAL_LINKED_WITH_TBB
# ifdef _DEBUG
  tbb::task_scheduler_init init(1);
# else
  tbb::task_scheduler_init init(10);
# endif
#endif

  int i = 0;
  bool stop = false;
  //for ( ; !stop ; ++i)
  {
    Wall_clock_timer t;
    CGAL::default_random = CGAL::Random(i);
    std::cerr << "Random seed = " << i << std::endl;
  
    //std::vector<Point> points = generate_points_on_plane<Point>(NUM_POINTS);
    //std::vector<Point> points = generate_points_on_sphere_3<Point>(NUM_POINTS, 3.0);
    std::vector<Point> points = generate_points_on_sphere_d<Point>(NUM_POINTS, AMBIENT_DIMENSION, 3.0);
    //std::vector<Point> points = generate_points_on_klein_bottle_3D<Point>(NUM_POINTS, 4., 3.);
    //std::vector<Point> points = generate_points_on_klein_bottle_4D<Point>(NUM_POINTS, 4., 3.);

    CGAL::Tangential_complex<
      Kernel, 
      INTRINSIC_DIMENSION, 
      CGAL::Parallel_tag> tc(points.begin(), points.end());
    double init_time = t.elapsed(); t.reset();

    tc.compute_tangential_complex();
    double computation_time = t.elapsed(); t.reset();

    std::set<std::set<std::size_t> > incorrect_simplices;
    //stop = !tc.check_if_all_simplices_are_in_the_ambient_delaunay(&incorrect_simplices);

    std::stringstream output_filename;
    output_filename << "data/test_tc_" << INTRINSIC_DIMENSION
      << "_in_R" << AMBIENT_DIMENSION << ".off";
    std::ofstream off_stream(output_filename.str());
    if (INTRINSIC_DIMENSION <= 3)
      tc.export_to_off(off_stream, true, &incorrect_simplices, true);
    else
      tc.number_of_inconsistent_simplices();

    double export_time = t.elapsed(); t.reset();

    std::cerr << std::endl
      << "================================================" << std::endl
      << "Computation times (seconds): " << std::endl
      << "  * Tangential complex: " << init_time + computation_time
      << std::endl
      << "    - Init + kd-tree = " << init_time << std::endl
      << "    - TC computation = " << computation_time << std::endl
      << "  * Export to OFF: " << export_time << std::endl
      << "================================================" << std::endl
      << std::endl;
  }

  return 0;
}