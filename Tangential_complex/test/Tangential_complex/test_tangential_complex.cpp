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
#ifdef CGAL_LINKED_WITH_TBB
# ifdef _DEBUG
  tbb::task_scheduler_init init(1);
# else
  tbb::task_scheduler_init init(10);
# endif
#endif

  const int INTRINSIC_DIMENSION = 1;
  const int AMBIENT_DIMENSION = 3;

  typedef CGAL::Epick_d<CGAL::Dimension_tag<AMBIENT_DIMENSION> > Kernel;
  typedef Kernel::Point_d                                        Point;
 
  int i = 0;
  bool stop = false;
  //for ( ; !stop ; ++i)
  {
    Wall_clock_timer t;
    CGAL::default_random = CGAL::Random(i);
    std::cerr << "Random seed = " << i << std::endl;
  
    std::vector<Point> points = 
      //generate_points_on_circle_2<Kernel>(NUM_POINTS, 3.);
      generate_points_on_moment_curve<Kernel>(NUM_POINTS, AMBIENT_DIMENSION, 0., 1.);
      //generate_points_on_plane<Kernel>(NUM_POINTS);
      //generate_points_on_sphere_3<Kernel>(NUM_POINTS, 3.0);
      //generate_points_on_sphere_d<Kernel>(NUM_POINTS, AMBIENT_DIMENSION, 3.0);
      //generate_points_on_klein_bottle_3D<Kernel>(NUM_POINTS, 4., 3.);
      //generate_points_on_klein_bottle_4D<Kernel>(NUM_POINTS, 4., 3.);

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