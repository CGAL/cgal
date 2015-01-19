//#undef CGAL_LINKED_WITH_TBB // CJTODO TEMP

// Without TBB_USE_THREADING_TOOL Intel Inspector XE will report false positives in Intel TBB
// (http://software.intel.com/en-us/articles/compiler-settings-for-threading-error-analysis-in-intel-inspector-xe/)
#ifdef _DEBUG
# define TBB_USE_THREADING_TOOL
#endif

#include <CGAL/Epick_d.h>
#include <CGAL/Tangential_complex.h>
#include <CGAL/Random.h>
#include <CGAL/Mesh_3/Profiling_tools.h>

#include "test_utilities.h"

#include <fstream>
#include <math.h>

#ifdef CGAL_LINKED_WITH_TBB
# include <tbb/task_scheduler_init.h>
#endif

#ifdef _DEBUG
  const int NUM_POINTS = 50;
#else
  const int NUM_POINTS = 30000;
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

  const int INTRINSIC_DIMENSION = 3;
  const int AMBIENT_DIMENSION   = 9;

  typedef CGAL::Epick_d<CGAL::Dimension_tag<AMBIENT_DIMENSION> > Kernel;
  typedef Kernel::FT                                             FT;
  typedef Kernel::Point_d                                        Point;
 
  int i = 0;
  bool stop = false;
  //for ( ; !stop ; ++i)
  {
    Kernel k;
    Wall_clock_timer t;
    CGAL::default_random = CGAL::Random(i);
    std::cerr << "Random seed = " << i << std::endl;
    
#ifdef CGAL_TC_PROFILING
    Wall_clock_timer t_gen;
#endif

    /*std::vector<Point> points =
      //generate_points_on_circle_2<Kernel>(NUM_POINTS, 3.);
      //generate_points_on_moment_curve<Kernel>(NUM_POINTS, AMBIENT_DIMENSION, 0., 1.);
      //generate_points_on_plane<Kernel>(NUM_POINTS);
      //generate_points_on_sphere_3<Kernel>(NUM_POINTS, 3.0);
      //generate_points_on_sphere_d<Kernel>(NUM_POINTS, AMBIENT_DIMENSION, 3.0);
      //generate_points_on_klein_bottle_3D<Kernel>(NUM_POINTS, 4., 3.);
      generate_points_on_klein_bottle_4D<Kernel>(NUM_POINTS, 4., 3.);
      //generate_points_on_klein_bottle_variant_5D<Kernel>(NUM_POINTS, 4., 3.);*/

    // LOAD FROM A FILE
    std::vector<Point> points;
    load_points_from_file<Point>(
      "data/SO3_10000.txt", std::back_inserter(points));

#ifdef CGAL_TC_PROFILING
    std::cerr << "Point set generated in " << t_gen.elapsed()
              << " seconds." << std::endl;
#endif

    std::size_t num_points_before = points.size();
    points = sparsify_point_set(
      k, points, FT(INPUT_SPARSITY)*FT(INPUT_SPARSITY));
    std::cerr << "Number of points before/after sparsification: "
      << num_points_before << " / " << points.size() << std::endl;

    CGAL::Tangential_complex<
      Kernel,
      INTRINSIC_DIMENSION,
      CGAL::Parallel_tag> tc(points.begin(), points.end(), k);
    double init_time = t.elapsed(); t.reset();

    tc.compute_tangential_complex();
    double computation_time = t.elapsed(); t.reset();
    
    std::set<std::set<std::size_t> > incorrect_simplices;
    //stop = !tc.check_if_all_simplices_are_in_the_ambient_delaunay(&incorrect_simplices);

    double export_before_time = -1.;
    if (INTRINSIC_DIMENSION <= 3)
    {
      t.reset();
      std::stringstream output_filename;
      output_filename << "output/test_tc_" << INTRINSIC_DIMENSION
        << "_in_R" << AMBIENT_DIMENSION << "_BEFORE_FIX.off";
      std::ofstream off_stream(output_filename.str().c_str());
      tc.export_to_off(off_stream, true, &incorrect_simplices, true);
      export_before_time = t.elapsed(); t.reset();
    }


    t.reset();
    tc.fix_inconsistencies();
    double fix_time = t.elapsed(); t.reset();

    double export_after_time = -1.;
    if (INTRINSIC_DIMENSION <= 3)
    {
      t.reset();
      std::stringstream output_filename;
      output_filename << "output/test_tc_" << INTRINSIC_DIMENSION
        << "_in_R" << AMBIENT_DIMENSION << "_AFTER_FIX.off";
      std::ofstream off_stream(output_filename.str().c_str());
      tc.export_to_off(off_stream, true, &incorrect_simplices, true);
      export_after_time = t.elapsed(); t.reset();
    }
    /*else
      tc.number_of_inconsistent_simplices();*/


    std::cerr << std::endl
      << "================================================" << std::endl
      << "Number of vertices: " << tc.number_of_vertices() << std::endl
      << "Computation times (seconds): " << std::endl
      << "  * Tangential complex: " << init_time + computation_time
      << std::endl
      << "    - Init + kd-tree = " << init_time << std::endl
      << "    - TC computation = " << computation_time << std::endl
      << "  * Export to OFF (before fix): " << export_before_time << std::endl
      << "  * Fix inconsistencies: " << fix_time << std::endl
      << "  * Export to OFF (after fix): " << export_after_time << std::endl
      << "================================================" << std::endl
      << std::endl;
  }

  return 0;
}
