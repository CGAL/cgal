//#undef CGAL_LINKED_WITH_TBB // CJTODO TEMP

// Without TBB_USE_THREADING_TOOL Intel Inspector XE will report false positives in Intel TBB
// (http://software.intel.com/en-us/articles/compiler-settings-for-threading-error-analysis-in-intel-inspector-xe/)
#ifdef _DEBUG
# define TBB_USE_THREADING_TOOL
#endif

#include <CGAL/assertions_behaviour.h>
#include <CGAL/Epick_d.h>
#include <CGAL/Tangential_complex.h>
#include <CGAL/Random.h>
#include <CGAL/Mesh_3/Profiling_tools.h>

#include "testing_utilities.h"

#include <fstream>
#include <math.h>

#ifdef CGAL_LINKED_WITH_TBB
# include <tbb/task_scheduler_init.h>
#endif

//=============== Constants =================
const double INPUT_SPARSITY = 0.05;
#ifdef _DEBUG
  const int NUM_POINTS = 50;
#else
  const int NUM_POINTS = 30000;
#endif
//===========================================

int main()
{
#if defined(CHECK_MEMORY_LEAKS_ON_MSVC) && defined(_MSC_VER)
  _CrtSetDbgFlag ( _CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF );
#endif
  CGAL::set_error_behaviour(CGAL::ABORT);

#ifdef CGAL_LINKED_WITH_TBB
# ifdef _DEBUG
  tbb::task_scheduler_init init(1);
# else
  tbb::task_scheduler_init init(10);
# endif
#endif

  const int INTRINSIC_DIMENSION = 3;
  const int AMBIENT_DIMENSION   = 9;

  typedef CGAL::Epick_d<CGAL::Dimension_tag<AMBIENT_DIMENSION> >  Kernel;
  typedef Kernel::FT                                              FT;
  typedef Kernel::Point_d                                         Point;
  typedef CGAL::Tangential_complex<
    Kernel, CGAL::Dimension_tag<INTRINSIC_DIMENSION>,
    CGAL::Parallel_tag>                                           TC;

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

    TC tc(points.begin(), points.end(), INPUT_SPARSITY, INTRINSIC_DIMENSION, k);
    double init_time = t.elapsed(); t.reset();

    tc.compute_tangential_complex();
    double computation_time = t.elapsed(); t.reset();
    
    if (ambient_dim <= 4)
      tc.check_if_all_simplices_are_in_the_ambient_delaunay();

    double export_before_time = -1.;
    if (INTRINSIC_DIMENSION <= 3)
    {
      t.reset();
      std::stringstream output_filename;
      output_filename << "output/test_tc_" << INTRINSIC_DIMENSION
        << "_in_R" << AMBIENT_DIMENSION << "_BEFORE_FIX.off";
      std::ofstream off_stream(output_filename.str().c_str());
      tc.export_to_off(off_stream, true);
      export_before_time = t.elapsed(); t.reset();
    }


    t.reset();
    unsigned int num_fix_steps;
    std::size_t initial_num_inconsistent_local_tr;
    std::size_t best_num_inconsistent_local_tr;
    std::size_t final_num_inconsistent_local_tr;
    CGAL::Fix_inconsistencies_status fix_ret = 
      tc.fix_inconsistencies_using_perturbation(
        num_fix_steps, initial_num_inconsistent_local_tr,
        best_num_inconsistent_local_tr, final_num_inconsistent_local_tr, 1000.);
    double fix_time = t.elapsed(); t.reset();

    double export_after_time = -1.;
    if (INTRINSIC_DIMENSION <= 3)
    {
      t.reset();
      std::stringstream output_filename;
      output_filename << "output/test_tc_" << INTRINSIC_DIMENSION
        << "_in_R" << AMBIENT_DIMENSION << "_AFTER_FIX.off";
      std::ofstream off_stream(output_filename.str().c_str());
      tc.export_to_off(off_stream, true);
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
      << "  * Fix inconsistencies: " << fix_time
      <<      " (" << num_fix_steps << " steps) ==> "
      <<      (fix_ret == CGAL::TC_FIXED ? "FIXED" : "NOT fixed") << std::endl
      << "  * Export to OFF (after fix): " << export_after_time << std::endl
      << "================================================" << std::endl
      << std::endl;
  }

  return 0;
}
