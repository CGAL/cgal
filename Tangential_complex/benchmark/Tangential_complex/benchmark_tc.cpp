//#undef CGAL_LINKED_WITH_TBB // CJTODO TEMP

// Without TBB_USE_THREADING_TOOL Intel Inspector XE will report false positives in Intel TBB
// (http://software.intel.com/en-us/articles/compiler-settings-for-threading-error-analysis-in-intel-inspector-xe/)
#ifdef _DEBUG
# define TBB_USE_THREADING_TOOL
#endif

#include <cstddef>

//#define CGAL_TC_USE_ANOTHER_POINT_SET_FOR_TANGENT_SPACE_ESTIM
//#define TC_PROTECT_POINT_SET_DELTA  0.003
//#define JUST_BENCHMARK_SPATIAL_SEARCH // CJTODO: test
//#define CHECK_IF_ALL_SIMPLICES_ARE_IN_THE_AMBIENT_DELAUNAY
//#define TC_INPUT_STRIDES 3 // only take one point every TC_INPUT_STRIDES points
//#define TC_NO_EXPORT
//#define CGAL_TC_ALVAREZ_SURFACE_WINDOW 0.95 // 1.9 - 0.95
//#define CGAL_TC_EXPORT_SPARSIFIED_POINT_SET
//#define CGAL_TC_EXPORT_ALL_COORDS_IN_OFF

const std::size_t ONLY_LOAD_THE_FIRST_N_POINTS = 500000; // 1e10

#include <CGAL/assertions_behaviour.h>
#include <CGAL/Epick_d.h>
#include <CGAL/Tangential_complex.h>
#include <CGAL/Random.h>
#include <CGAL/Mesh_3/Profiling_tools.h>

#include "../../test/Tangential_complex/testing_utilities.h"

#include <boost/algorithm/string/replace.hpp>
#include <boost/algorithm/string/trim_all.hpp>
#include <boost/range/adaptor/strided.hpp>

#include <cstdlib>
#include <ctime>
#include <fstream>
#include <math.h>

#ifdef CGAL_LINKED_WITH_TBB
# include <tbb/task_scheduler_init.h>
#endif
#include "XML_exporter.h"
#define CGAL_TC_EXPORT_PERFORMANCE_DATA
#define CGAL_TC_SET_PERFORMANCE_DATA(value_name, value) \
        XML_perf_data::set(value_name, value);

const char * const BENCHMARK_SCRIPT_FILENAME = "benchmark_script.txt";

typedef CGAL::Epick_d<CGAL::Dynamic_dimension_tag>              Kernel;
typedef Kernel::FT                                              FT;
typedef Kernel::Point_d                                         Point;
typedef Kernel::Vector_d                                        Vector;
typedef CGAL::Tangential_complex<
  Kernel, CGAL::Dynamic_dimension_tag,
  CGAL::Parallel_tag>                                           TC;


#ifdef TC_PROTECT_POINT_SET_DELTA
# include <CGAL/Tangential_complex/protected_sets.h> // CJTODO TEST
#endif

#ifdef JUST_BENCHMARK_SPATIAL_SEARCH
std::ofstream spatial_search_csv_file("benchmark_spatial_search.csv");
#endif

using namespace CGAL::Tangential_complex_;

class XML_perf_data
{
public:
  typedef Streaming_XML_exporter<std::string> XML_exporter;

  XML_perf_data(const std::string &filename)
    : m_xml(filename, "ContainerPerformance", "Perf",
            construct_subelements_names())
  {}

  virtual ~XML_perf_data()
  {
  }

  static XML_perf_data &get()
  {
    static XML_perf_data singleton(build_filename());
    return singleton;
  }

  template <typename Value_type>
  static void set(const std::string &name, Value_type value)
  {
    get().set_data(name, value);
  }

  static void commit()
  {
    get().commit_current_element();
  }

protected:
  static std::string build_filename()
  {
    std::stringstream sstr;
    sstr << "perf_logs/Performance_log_" << time(0) << ".xml";
    return sstr.str();
  }

  static std::vector<std::string> construct_subelements_names()
  {
    std::vector<std::string> subelements;
    subelements.push_back("Input");
    subelements.push_back("Intrinsic_dim");
    subelements.push_back("Ambient_dim");
    subelements.push_back("Sparsity");
    subelements.push_back("Num_points_in_input");
    subelements.push_back("Num_points");
    subelements.push_back("Initial_num_inconsistent_local_tr");
    subelements.push_back("Best_num_inconsistent_local_tr");
    subelements.push_back("Final_num_inconsistent_local_tr");
    subelements.push_back("Init_time");
    subelements.push_back("Comput_time");
    subelements.push_back("Perturb_successful");
    subelements.push_back("Perturb_time");
    subelements.push_back("Perturb_steps");
    subelements.push_back("Add_higher_dim_simpl_time");
    subelements.push_back("Result_pure_pseudomanifold");
    subelements.push_back("Result_num_wrong_dim_simplices");
    subelements.push_back("Result_num_wrong_number_of_cofaces");
    subelements.push_back("Result_num_unconnected_stars");
    subelements.push_back("Info");

    return subelements;
  }

  void set_data(const std::string &name, const std::string &value)
  {
    m_current_element[name] = value;
  }

  template <typename Value_type>
  void set_data(const std::string &name, Value_type value)
  {
    std::stringstream sstr;
    sstr << value;
    set_data(name, sstr.str());
  }

  void commit_current_element()
  {
    m_xml.add_element(m_current_element);
    m_current_element.clear();
  }

  XML_exporter m_xml;
  XML_exporter::Element_with_map m_current_element;
};

class Test_dim
{
public:
  Test_dim(
    int min_allowed_dim = 0, 
    int max_allowed_dim = std::numeric_limits<int>::max())
    : m_min_allowed_dim(min_allowed_dim), m_max_allowed_dim(max_allowed_dim)
  {}

  template <typename Simplex>
  bool operator()(Simplex const& s)
  {
    return s.size() - 1 >= m_min_allowed_dim
      && s.size() - 1 <= m_max_allowed_dim;
  }

private:
  int m_min_allowed_dim;
  int m_max_allowed_dim;
};

// color_inconsistencies: only works if p_complex = NULL
template <typename TC>
bool export_to_off(
  TC const& tc, 
  std::string const& input_name_stripped,
  std::string const& suffix,
  bool color_inconsistencies = false,
  typename TC::Simplicial_complex const* p_complex = NULL,
  std::set<std::set<std::size_t> > const *p_simpl_to_color_in_red = NULL,
  std::set<std::set<std::size_t> > const *p_simpl_to_color_in_green = NULL,
  std::set<std::set<std::size_t> > const *p_simpl_to_color_in_blue = NULL)
{
#ifdef TC_NO_EXPORT
  return true;
#endif

#if 0
  Kernel k;
  FT center_pt[] = { -0.5, -CGAL::sqrt(3.) / 2, -0.5, CGAL::sqrt(3.) / 2 };
  FT proj_pt[] = { 0., 0., 0., 0.2 };
  S3_to_R3_stereographic_projection<Kernel>
    proj_functor(0.2, 
                 Point(4, &center_pt[0], &center_pt[4]),
                 k);
#else
  CGAL::Identity<Point> proj_functor;
  //Kernel k;
  //std::array<int, 3> sel = { 1, 3, 5 };
  //Orthogonal_projection<Kernel> proj_functor(sel, k);
#endif

  if (tc.intrinsic_dimension() <= 3)
  {
    std::stringstream output_filename;
    output_filename << "output/" << input_name_stripped << "_" 
      << tc.intrinsic_dimension() << "_in_R" 
      << tc.ambient_dimension() << "_"
      << tc.number_of_vertices() << "v"
      << suffix << ".off";
    std::ofstream off_stream(output_filename.str().c_str());

    if (p_complex)
    {
#ifndef TC_NO_EXPORT
      tc.export_to_off(
        *p_complex, off_stream, 
        p_simpl_to_color_in_red,
        p_simpl_to_color_in_green, 
        p_simpl_to_color_in_blue,
        proj_functor);
#endif
    }
    else
    {
/*#ifdef CGAL_ALPHA_TC
      TC::Simplicial_complex complex;
      tc.export_TC(complex, false);
      tc.export_to_off(
        complex, off_stream, 
        p_simpl_to_color_in_red,
        p_simpl_to_color_in_green, 
        p_simpl_to_color_in_blue);
#else*/
      tc.export_to_off(
        off_stream, color_inconsistencies, 
        p_simpl_to_color_in_red,
        p_simpl_to_color_in_green, 
        p_simpl_to_color_in_blue,
        NULL,
        proj_functor);
//#endif
    }
    return true;
  }
  return false;
}

void make_tc(std::vector<Point> &points, 
             TC::TS_container const& tangent_spaces, // can be empty
             int intrinsic_dim,
             bool sparsify = true,
             double sparsity = 0.01, 
             bool perturb = true, 
             bool add_high_dim_simpl = false, 
             bool collapse = false,
             double time_limit_for_perturb = 0.,
             const char *input_name = "tc")
{
  Kernel k;
  
  if (sparsify && !tangent_spaces.empty())
  {
    std::cerr << "Error: cannot sparsify point set with pre-computed normals.\n";
    return;
  }

  // CJTODO TEMP TEST
  //TC::Simplicial_complex compl;
  //{std::size_t ss[] = {0, 1, 2}; compl.add_simplex(std::set<std::size_t>(ss, ss + 3)); }
  //{std::size_t ss[] = {0, 2, 3}; compl.add_simplex(std::set<std::size_t>(ss, ss + 3)); }
  //{std::size_t ss[] = {0, 3, 4}; compl.add_simplex(std::set<std::size_t>(ss, ss + 3)); }
  //{std::size_t ss[] = {0, 4, 1}; compl.add_simplex(std::set<std::size_t>(ss, ss + 3)); }
  //{std::size_t ss[] = {0, 5, 6}; compl.add_simplex(std::set<std::size_t>(ss, ss + 3)); }
  //compl.is_pure_pseudomanifold(2, 7, false, 10);

  //TC::Simplicial_complex compl;
  //{std::size_t ss[] = {0, 1, 2, 5}; compl.add_simplex(std::set<std::size_t>(ss, ss + 4)); }
  //{std::size_t ss[] = {0, 2, 3, 5}; compl.add_simplex(std::set<std::size_t>(ss, ss + 4)); }
  //{std::size_t ss[] = {0, 3, 4, 5}; compl.add_simplex(std::set<std::size_t>(ss, ss + 4)); }
  //{std::size_t ss[] = {0, 4, 1, 5}; compl.add_simplex(std::set<std::size_t>(ss, ss + 4)); }
  //{std::size_t ss[] = {0, 1, 2, 6}; compl.add_simplex(std::set<std::size_t>(ss, ss + 4)); }
  //{std::size_t ss[] = {0, 2, 3, 6}; compl.add_simplex(std::set<std::size_t>(ss, ss + 4)); }
  //{std::size_t ss[] = {0, 3, 4, 6}; compl.add_simplex(std::set<std::size_t>(ss, ss + 4)); }
  //{std::size_t ss[] = {0, 4, 1, 6}; compl.add_simplex(std::set<std::size_t>(ss, ss + 4)); }
  //{std::size_t ss[] = {0, 4, 7, 8}; compl.add_simplex(std::set<std::size_t>(ss, ss + 4)); }
  //compl.is_pure_pseudomanifold(3, 9, false, 10);
  // /CJTODO TEMP TEST

#ifdef JUST_BENCHMARK_SPATIAL_SEARCH
  benchmark_spatial_search(points, k, spatial_search_csv_file);
  return;
#endif

  //===========================================================================
  // Init
  //===========================================================================
  Wall_clock_timer t;

  // Get input_name_stripped
  std::string input_name_stripped(input_name);
  size_t slash_index = input_name_stripped.find_last_of('/');
  if (slash_index == std::string::npos)
    slash_index = input_name_stripped.find_last_of('\\');
  if (slash_index == std::string::npos)
    slash_index = 0;
  else
    ++slash_index;
  input_name_stripped = input_name_stripped.substr(
    slash_index, input_name_stripped.find_last_of('.') - slash_index);

  int ambient_dim = k.point_dimension_d_object()(*points.begin());

  CGAL_TC_SET_PERFORMANCE_DATA("Num_points_in_input", points.size());

#ifdef CGAL_TC_USE_ANOTHER_POINT_SET_FOR_TANGENT_SPACE_ESTIM
  std::vector<Point> points_not_sparse = points;
#endif

  //===========================================================================
  // Sparsify point set if requested
  //===========================================================================
  if (sparsify)
  {
    std::size_t num_points_before = points.size();
    points = sparsify_point_set(k, points, sparsity*sparsity);
    std::cerr << "Number of points before/after sparsification: "
      << num_points_before << " / " << points.size() << "\n";

#ifdef CGAL_TC_EXPORT_SPARSIFIED_POINT_SET
    std::ofstream ps_stream("output/sparsified_point_set.txt");
    export_point_set(k, points, ps_stream);
#endif
  }

#ifdef TC_PROTECT_POINT_SET_DELTA
  // CJTODO TEST    
# ifdef CGAL_TC_PROFILING
  Wall_clock_timer t_protection;
# endif

  std::vector<Point> points2;
  std::vector<int> dummy;
  std::vector<std::vector<int> > dummy2;
  landmark_choice_protected_delaunay(
    points, points.size(), points2, dummy, TC_PROTECT_POINT_SET_DELTA, dummy2, false, true);
  points = points2;

# ifdef CGAL_TC_PROFILING
  std::cerr << "Point set protected in " << t_protection.elapsed()
    << " seconds.\n";
# endif

  std::cerr << "Number of points after PROTECTION: " << points.size() << "\n";
#endif

  CGAL_TC_SET_PERFORMANCE_DATA("Sparsity", sparsity);
  CGAL_TC_SET_PERFORMANCE_DATA("Num_points", points.size());

  //===========================================================================
  // Compute Tangential Complex
  //===========================================================================

#ifdef CGAL_TC_USE_ANOTHER_POINT_SET_FOR_TANGENT_SPACE_ESTIM
  TC tc(points.begin(), points.end(), sparsity, intrinsic_dim,
    points_not_sparse.begin(), points_not_sparse.end(), k);
#else
  TC tc(points.begin(), points.end(), sparsity, intrinsic_dim, k);
#endif

  if (!tangent_spaces.empty())
  {
    tc.set_tangent_planes(tangent_spaces);
  }

  double init_time = t.elapsed(); t.reset();

  tc.compute_tangential_complex();
  double computation_time = t.elapsed(); t.reset();

#ifdef CHECK_IF_ALL_SIMPLICES_ARE_IN_THE_AMBIENT_DELAUNAY
  if (ambient_dim <= 4)
    tc.check_if_all_simplices_are_in_the_ambient_delaunay();
#endif

  //tc.check_correlation_between_inconsistencies_and_fatness();

  //===========================================================================
  // Export to OFF
  //===========================================================================

  // Create complex
  int max_dim = -1;
  TC::Simplicial_complex complex;
  std::set<std::set<std::size_t> > inconsistent_simplices;
  max_dim = tc.export_TC(complex, false, 2, &inconsistent_simplices);

  t.reset();
  double export_before_time = 
    (export_to_off(tc, input_name_stripped, "_INITIAL_TC", true, 
      &complex, &inconsistent_simplices) ? t.elapsed() : -1);
  t.reset();
  
  unsigned int num_perturb_steps = 0;
  double perturb_time = -1;
    double export_after_perturb_time = -1.;
  CGAL::Fix_inconsistencies_status perturb_ret = CGAL::FIX_NOT_PERFORMED;
  if (perturb)
  {
    //=========================================================================
    // Try to fix inconsistencies by perturbing points
    //=========================================================================
    t.reset();
    std::size_t initial_num_inconsistent_local_tr;
    std::size_t best_num_inconsistent_local_tr;
    std::size_t final_num_inconsistent_local_tr;
    perturb_ret = tc.fix_inconsistencies_using_perturbation(
      num_perturb_steps, initial_num_inconsistent_local_tr,
      best_num_inconsistent_local_tr, final_num_inconsistent_local_tr,
      time_limit_for_perturb);
    perturb_time = t.elapsed(); t.reset();

    CGAL_TC_SET_PERFORMANCE_DATA("Initial_num_inconsistent_local_tr", 
                                 initial_num_inconsistent_local_tr);
    CGAL_TC_SET_PERFORMANCE_DATA("Best_num_inconsistent_local_tr", 
                                 best_num_inconsistent_local_tr);
    CGAL_TC_SET_PERFORMANCE_DATA("Final_num_inconsistent_local_tr", 
                                 final_num_inconsistent_local_tr);

    //tc.check_correlation_between_inconsistencies_and_fatness();

    // DEBUGGING: confirm that all stars were actually refreshed
    //std::cerr << yellow << "FINAL CHECK...\n" << white;
    //std::size_t num_inc = tc.number_of_inconsistent_simplices(true).second;
    //tc.refresh_tangential_complex();
    //if (CGAL::cpp11::get<1>(tc.number_of_inconsistent_simplices(true)) != num_inc)
    //  std::cerr << red << "FINAL CHECK: FAILED.\n" << white;
    //else
    //  std::cerr << green << "FINAL CHECK: PASSED.\n" << white;


    //=========================================================================
    // Export to OFF
    //=========================================================================

    // Re-build the complex
    std::set<std::set<std::size_t> > inconsistent_simplices;
    max_dim = tc.export_TC(complex, false, 2, &inconsistent_simplices);

    t.reset();
    bool exported = export_to_off(
      tc, input_name_stripped, "_AFTER_FIX", true, &complex, 
      &inconsistent_simplices);
    double export_after_perturb_time = (exported ? t.elapsed() : -1);
    t.reset();

    //std::string fn = "output/inc_stars/";
    //fn += input_name_stripped;
    //tc.export_inconsistent_stars_to_OFF_files(fn);
  }
  else
  {
    CGAL_TC_SET_PERFORMANCE_DATA("Initial_num_inconsistent_local_tr", "N/A");
    CGAL_TC_SET_PERFORMANCE_DATA("Best_num_inconsistent_local_tr", "N/A");
    CGAL_TC_SET_PERFORMANCE_DATA("Final_num_inconsistent_local_tr", "N/A");
  }

  // CJTODO TEST
  //tc.check_and_solve_inconsistencies_by_filtering_simplices_out();

  double fix2_time = -1;
  double export_after_fix2_time = -1.;
  if (add_high_dim_simpl)
  {
    //=========================================================================
    // Try to fix inconsistencies by adding higher-dimension simplices
    //=========================================================================
    t.reset();
    // Try to solve the remaining inconstencies
#ifdef CGAL_ALPHA_TC
    tc.solve_inconsistencies_using_alpha_TC();
#else
    tc.check_and_solve_inconsistencies_by_adding_higher_dim_simplices();
#endif
    fix2_time = t.elapsed(); t.reset();

    /*std::set<std::set<std::size_t> > not_delaunay_simplices;
    if (ambient_dim <= 4)
    {
      tc.check_if_all_simplices_are_in_the_ambient_delaunay(
        &complex, true, &not_delaunay_simplices);
    }*/
  
    //=========================================================================
    // Export to OFF
    //=========================================================================

    // Re-build the complex
    std::set<std::set<std::size_t> > inconsistent_simplices;
    max_dim = tc.export_TC(complex, false, 2, &inconsistent_simplices);

    t.reset();
    bool exported = export_to_off(
      tc, input_name_stripped, "_AFTER_FIX2", false, &complex, 
      &inconsistent_simplices);
    double export_after_fix2_time = (exported ? t.elapsed() : -1);
    t.reset();
  }
  else
  {
    std::set<std::set<std::size_t> > inconsistent_simplices;
    max_dim = tc.export_TC(complex, false, 2, &inconsistent_simplices);
  }

  complex.display_stats();

  if (intrinsic_dim == 2)
    complex.euler_characteristic(true);

  // CJTODO TEMP: Export to OFF with higher-dim simplices colored
  /*std::set<std::set<std::size_t> > higher_dim_simplices;
  complex.get_simplices_matching_test(
    Test_dim(intrinsic_dim + 1),
    std::inserter(higher_dim_simplices, higher_dim_simplices.begin()));
  export_to_off(
    tc, input_name_stripped, "_BEFORE_COLLAPSE", false, &complex, 
    &inconsistent_simplices, &higher_dim_simplices);*/
  
  //===========================================================================
  // Collapse
  //===========================================================================
  if (collapse)
  {
    complex.collapse(max_dim);
    complex.display_stats();
  }

  //===========================================================================
  // Is the result a pure pseudomanifold?
  //===========================================================================
  std::size_t num_wrong_dim_simplices, 
              num_wrong_number_of_cofaces, 
              num_unconnected_stars;
  std::set<std::set<std::size_t> > wrong_dim_simplices;
  std::set<std::set<std::size_t> > wrong_number_of_cofaces_simplices;
  std::set<std::set<std::size_t> > unconnected_stars_simplices;
  bool is_pure_pseudomanifold = complex.is_pure_pseudomanifold(
    intrinsic_dim, tc.number_of_vertices(), 
#ifdef CGAL_TC_ALVAREZ_SURFACE_WINDOW
    true, // allow borders
#else
    false, // do NOT allow borders
#endif
    false, 1,
    &num_wrong_dim_simplices, &num_wrong_number_of_cofaces, 
    &num_unconnected_stars,
    &wrong_dim_simplices, &wrong_number_of_cofaces_simplices, 
    &unconnected_stars_simplices);

  //===========================================================================
  // Export to OFF
  //===========================================================================
  
  double export_after_collapse_time = -1.;
  if (collapse)
  {
    t.reset();
    bool exported = export_to_off(
      tc, input_name_stripped, "_AFTER_COLLAPSE", false, &complex,
      &wrong_dim_simplices, &wrong_number_of_cofaces_simplices,
      &unconnected_stars_simplices);
    std::cerr
      << " OFF colors:\n"
      << "   * Red: wrong dim simplices\n"
      << "   * Green: wrong number of cofaces simplices\n"
      << "   * Blue: not-connected stars\n";
    double export_after_collapse_time = (exported ? t.elapsed() : -1.);
    t.reset();
  }

  //===========================================================================
  // Display info
  //===========================================================================

  std::cerr
    << "\n================================================\n"
    << "Number of vertices: " << tc.number_of_vertices() << "\n"
    << "Computation times (seconds): \n"
    << "  * Tangential complex: " << init_time + computation_time << "\n"
    << "    - Init + kd-tree = " << init_time << "\n"
    << "    - TC computation = " << computation_time << "\n"
    << "  * Export to OFF (before perturb): " << export_before_time << "\n"
    << "  * Fix inconsistencies 1: " << perturb_time
    <<      " (" << num_perturb_steps << " steps) ==> "
    <<      (perturb_ret == CGAL::TC_FIXED ? "FIXED" : "NOT fixed") << "\n"
    << "  * Fix inconsistencies 2: " << fix2_time << "\n"
    << "  * Export to OFF (after perturb): " << export_after_perturb_time << "\n"
    << "  * Export to OFF (after fix2): "<< export_after_fix2_time << "\n"
    << "  * Export to OFF (after collapse): "
    <<      export_after_collapse_time << "\n"
    << "================================================\n";
  
  //===========================================================================
  // Export info
  //===========================================================================
  CGAL_TC_SET_PERFORMANCE_DATA("Init_time", init_time);
  CGAL_TC_SET_PERFORMANCE_DATA("Comput_time", computation_time);
  CGAL_TC_SET_PERFORMANCE_DATA("Perturb_successful",
                                (perturb_ret == CGAL::TC_FIXED ? "Y" : "N"));
  CGAL_TC_SET_PERFORMANCE_DATA("Perturb_time", perturb_time);
  CGAL_TC_SET_PERFORMANCE_DATA("Perturb_steps", num_perturb_steps);
  CGAL_TC_SET_PERFORMANCE_DATA("Add_higher_dim_simpl_time", fix2_time);
  CGAL_TC_SET_PERFORMANCE_DATA("Result_pure_pseudomanifold",
                                (is_pure_pseudomanifold ? "Y" : "N"));
  CGAL_TC_SET_PERFORMANCE_DATA("Result_num_wrong_dim_simplices",
                                num_wrong_dim_simplices);
  CGAL_TC_SET_PERFORMANCE_DATA("Result_num_wrong_number_of_cofaces", 
                                num_wrong_number_of_cofaces);
  CGAL_TC_SET_PERFORMANCE_DATA("Result_num_unconnected_stars", 
                                num_unconnected_stars);
  CGAL_TC_SET_PERFORMANCE_DATA("Info", "");
}

int main()
{
  CGAL::set_error_behaviour(CGAL::ABORT);

#ifdef CGAL_LINKED_WITH_TBB
# ifdef _DEBUG
  int num_threads = 1;
# else
  int num_threads = 8;
# endif
#endif

  unsigned int seed = static_cast<unsigned int>(time(NULL));
  CGAL::default_random = CGAL::Random(seed);
  std::cerr << "Random seed = " << seed << "\n";

  std::ifstream script_file;
  script_file.open(BENCHMARK_SCRIPT_FILENAME);
  // Script?
  // Script file format: each line gives
  //    - Filename (point set) or "generate_XXX" (point set generation)
  //    - Ambient dim
  //    - Intrinsic dim
  //    - Number of iterations with these parameters
  if (script_file.is_open())
  {
    int i = 1;
#ifdef CGAL_LINKED_WITH_TBB
# ifdef BENCHMARK_WITH_1_TO_MAX_THREADS
    for(num_threads = 1 ;
          num_threads <= tbb::task_scheduler_init::default_num_threads() ;
          ++num_threads)
# endif
#endif
    /*for (Concurrent_mesher_config::get().num_work_items_per_batch = 5 ;
      Concurrent_mesher_config::get().num_work_items_per_batch < 100 ;
      Concurrent_mesher_config::get().num_work_items_per_batch += 5)*/
    {
#ifdef CGAL_LINKED_WITH_TBB
      tbb::task_scheduler_init init(
        num_threads > 0 ? num_threads : tbb::task_scheduler_init::automatic);
#endif

      std::cerr << "Script file '" << BENCHMARK_SCRIPT_FILENAME << "' found.\n";
      script_file.seekg(0);
      while (script_file.good())
      {
        std::string line;
        std::getline(script_file, line);
        if (line.size() > 1 && line[0] != '#')
        {
          boost::replace_all(line, "\t", " ");
          boost::trim_all(line);
          std::cerr << "\n\n";
          std::cerr << "*****************************************\n";
          std::cerr << "******* " << line << "\n";
          std::cerr << "*****************************************\n";
          std::stringstream sstr(line);

          std::string input;
          std::string param1;
          std::string param2;
          std::string param3;
          std::size_t num_points;
          int ambient_dim;
          int intrinsic_dim;
          char sparsify;
          double sparsity;
          char perturb, add_high_dim_simpl, collapse;
          double time_limit_for_perturb;
          int num_iteration;
          sstr >> input;
          sstr >> param1;
          sstr >> param2;
          sstr >> param3;
          sstr >> num_points;
          sstr >> ambient_dim;
          sstr >> intrinsic_dim;
          sstr >> sparsify;
          sstr >> sparsity;
          sstr >> perturb;
          sstr >> add_high_dim_simpl;
          sstr >> collapse;
          sstr >> time_limit_for_perturb;
          sstr >> num_iteration;

          for (int j = 0 ; j < num_iteration ; ++j)
          {
            std::string input_stripped = input;
            size_t slash_index = input_stripped.find_last_of('/');
            if (slash_index == std::string::npos)
              slash_index = input_stripped.find_last_of('\\');
            if (slash_index == std::string::npos)
              slash_index = 0;
            else
              ++slash_index;
            input_stripped = input_stripped.substr(
              slash_index, input_stripped.find_last_of('.') - slash_index);

            CGAL_TC_SET_PERFORMANCE_DATA("Input", input_stripped);
            CGAL_TC_SET_PERFORMANCE_DATA("Ambient_dim", ambient_dim);
            CGAL_TC_SET_PERFORMANCE_DATA("Intrinsic_dim", intrinsic_dim);
#ifdef CGAL_LINKED_WITH_TBB
            CGAL_TC_SET_PERFORMANCE_DATA(
              "Num_threads",
              (num_threads == -1 ? tbb::task_scheduler_init::default_num_threads() : num_threads));
#else
            CGAL_TC_SET_PERFORMANCE_DATA("Num_threads", "N/A");
#endif

            std::cerr << "\nTC #" << i << "...\n";
          
#ifdef CGAL_TC_PROFILING
            Wall_clock_timer t_gen;
#endif

            std::vector<Point> points;
            TC::TS_container tangent_spaces;

            if (input == "generate_moment_curve")
            {
              points = generate_points_on_moment_curve<Kernel>(
                num_points, ambient_dim,
                std::atof(param1.c_str()), std::atof(param2.c_str()));
            }
            else if (input == "generate_plane")
            {
              points = generate_points_on_plane<Kernel>(
                num_points, intrinsic_dim, ambient_dim);
            }
            else if (input == "generate_sphere_d")
            {
              points = generate_points_on_sphere_d<Kernel>(
                num_points, ambient_dim,
                std::atof(param1.c_str()),  // radius
                std::atof(param2.c_str())); // radius_noise_percentage
            }
            else if (input == "generate_two_spheres_d")
            {
              points = generate_points_on_two_spheres_d<Kernel>(
                num_points, ambient_dim,
                std::atof(param1.c_str()),
                std::atof(param2.c_str()),
                std::atof(param3.c_str()));
            }
            else if (input == "generate_3sphere_and_circle_d")
            {
              CGAL_assertion(intrinsic_dim == 3);
              CGAL_assertion(ambient_dim == 5);
              points = generate_points_on_3sphere_and_circle<Kernel>(
                num_points,
                std::atof(param1.c_str()));
            }
            else if (input == "generate_torus_3D")
            {
              points = generate_points_on_torus_3D<Kernel>(
                num_points,
                std::atof(param1.c_str()),
                std::atof(param2.c_str()),
                param3 == "Y");
            }
            else if (input == "generate_torus_d")
            {
              points = generate_points_on_torus_d<Kernel>(
                num_points, 
                intrinsic_dim,
                param1 == "Y", // uniform
                std::atof(param2.c_str())); // radius_noise_percentage
            }
            else if (input == "generate_klein_bottle_3D")
            {
              points = generate_points_on_klein_bottle_3D<Kernel>(
                num_points,
                std::atof(param1.c_str()), std::atof(param2.c_str()));
            }
            else if (input == "generate_klein_bottle_4D")
            {
              points = generate_points_on_klein_bottle_4D<Kernel>(
                num_points,
                std::atof(param1.c_str()), std::atof(param2.c_str()));
            }
            else if (input == "generate_klein_bottle_variant_5D")
            {
              points = generate_points_on_klein_bottle_variant_5D<Kernel>(
                num_points,
                std::atof(param1.c_str()), std::atof(param2.c_str()));
            }
            else
            {
              load_points_from_file<Kernel, typename TC::Tangent_space_basis>(
                input, std::back_inserter(points), 
                std::back_inserter(tangent_spaces), 
                ONLY_LOAD_THE_FIRST_N_POINTS);
            }

#ifdef CGAL_TC_PROFILING
            std::cerr << "Point set generated/loaded in " << t_gen.elapsed()
                      << " seconds.\n";
#endif

            if (!points.empty())
            {
#if defined(TC_INPUT_STRIDES) && TC_INPUT_STRIDES > 1
              auto p = points | boost::adaptors::strided(TC_INPUT_STRIDES); // CJTODO C++11 (auto)
              std::vector<Point> points(p.begin(), p.end());
              std::cerr << "****************************************\n"
                << "WARNING: taking 1 point every " << TC_INPUT_STRIDES
                << " points.\n"
                << "****************************************\n";
#endif
              make_tc(points, tangent_spaces, intrinsic_dim, sparsify=='Y',
                      sparsity, perturb=='Y', add_high_dim_simpl=='Y', 
                      collapse=='Y', time_limit_for_perturb, input.c_str());

              std::cerr << "TC #" << i++ << " done.\n";
              std::cerr << "\n---------------------------------\n";
            }
            else
            {
              std::cerr << "TC #" << i++ << ": no points loaded.\n";
            }

            XML_perf_data::commit();
          }
        }
      }
      script_file.seekg(0);
      script_file.clear();
    }

    script_file.close();
  }
  // Or not script?
  else
  {
    std::cerr << "Script file '" << BENCHMARK_SCRIPT_FILENAME << "' NOT found.\n";
  }

  system("pause");
  return 0;
}
