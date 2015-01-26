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

#include "../../test/Tangential_complex/test_utilities.h"


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
typedef CGAL::Tangential_complex<
  Kernel, CGAL::Dynamic_dimension_tag, 
  CGAL::Parallel_tag>                                           TC;

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
    subelements.push_back("Num_points");
    subelements.push_back("Sparsity");
    subelements.push_back("Init_time");
    subelements.push_back("Comput_time");
    subelements.push_back("Fix_successful");
    subelements.push_back("Fix_time");
    subelements.push_back("Fix_steps");
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

const double NOISE = 0.01; 
#ifdef _DEBUG
  const int NUM_POINTS = 50;
#else
  const int NUM_POINTS = 30000;
#endif

void make_tc(std::vector<Point> &points, int intrinsic_dim,
             double sparsity = 0., double time_limit_for_fix = 0.)
{
  Kernel k;
  Wall_clock_timer t;

  int ambient_dim = k.point_dimension_d_object()(*points.begin());
    
#ifdef CGAL_TC_PROFILING
  Wall_clock_timer t_gen;
#endif
  
#ifdef CGAL_TC_PROFILING
  std::cerr << "Point set generated in " << t_gen.elapsed()
            << " seconds." << std::endl;
#endif

  if (sparsity != 0.)
  {
    std::size_t num_points_before = points.size();
    points = sparsify_point_set(
      k, points, FT(INPUT_SPARSITY)*FT(INPUT_SPARSITY));
    std::cerr << "Number of points before/after sparsification: "
      << num_points_before << " / " << points.size() << std::endl;
  }
  
  CGAL_TC_SET_PERFORMANCE_DATA("Num_points", points.size());
  CGAL_TC_SET_PERFORMANCE_DATA("Sparsity", sparsity);

  TC tc(points.begin(), points.end(), intrinsic_dim, k);
  double init_time = t.elapsed(); t.reset();

  tc.compute_tangential_complex();
  double computation_time = t.elapsed(); t.reset();
    
  std::set<std::set<std::size_t> > incorrect_simplices;
  //stop = !tc.check_if_all_simplices_are_in_the_ambient_delaunay(&incorrect_simplices);

  double export_before_time = -1.;
  if (intrinsic_dim <= 3)
  {
    t.reset();
    std::stringstream output_filename;
    output_filename << "output/test_tc_" << intrinsic_dim
      << "_in_R" << ambient_dim << "_BEFORE_FIX.off";
    std::ofstream off_stream(output_filename.str().c_str());
    tc.export_to_off(off_stream, true, &incorrect_simplices, true);
    export_before_time = t.elapsed(); t.reset();
  }


  t.reset();
  unsigned int num_fix_steps;
  CGAL::Fix_inconsistencies_status fix_ret =
    tc.fix_inconsistencies(num_fix_steps, time_limit_for_fix);
  double fix_time = t.elapsed(); t.reset();

  double export_after_time = -1.;
  if (intrinsic_dim <= 3)
  {
    t.reset();
    std::stringstream output_filename;
    output_filename << "output/test_tc_" << intrinsic_dim
      << "_in_R" << ambient_dim << "_AFTER_FIX.off";
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
    << "  * Fix inconsistencies: " << fix_time 
    <<      " (" << num_fix_steps << " steps) ==> " 
    <<      (fix_ret == CGAL::TC_FIXED ? "FIXED" : "NOT fixed") << std::endl
    << "  * Export to OFF (after fix): " << export_after_time << std::endl
    << "================================================" << std::endl
    << std::endl;

    CGAL_TC_SET_PERFORMANCE_DATA("Init_time",     init_time);
    CGAL_TC_SET_PERFORMANCE_DATA("Comput_time",   computation_time);
    CGAL_TC_SET_PERFORMANCE_DATA("Fix_successful",
                                 (fix_ret == CGAL::TC_FIXED ? "Y" : "N"));
    CGAL_TC_SET_PERFORMANCE_DATA("Fix_time",      fix_time);
    CGAL_TC_SET_PERFORMANCE_DATA("Fix_steps",     num_fix_steps);
    CGAL_TC_SET_PERFORMANCE_DATA("Info", "");
}

int main()
{
#ifdef CGAL_LINKED_WITH_TBB
# ifdef _DEBUG
  int num_threads = 1;
# else
  int num_threads = 10;
# endif
#endif

  int seed = CGAL::default_random.get_int(0, 1<<30);
  CGAL::default_random = CGAL::Random();
  std::cerr << "Random seed = " << seed << std::endl;
  
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

      std::cerr << "Script file '" << BENCHMARK_SCRIPT_FILENAME << "' found." << std::endl;
      script_file.seekg(0);
      while (script_file.good())
      {
        std::string line;
        std::getline(script_file, line);
        if (line.size() > 1 && line[0] != '#')
        {
          std::cerr << std::endl << std::endl;
          std::cerr << "*****************************************" << std::endl;
          std::cerr << "******* " << line << std::endl;
          std::cerr << "*****************************************" << std::endl;
          std::stringstream sstr(line);

          std::string input;
          int ambient_dim;
          int intrinsic_dim;
          double sparsity;
          double time_limit_for_fix;
          int num_iteration;
          sstr >> input;
          sstr >> ambient_dim;
          sstr >> intrinsic_dim;
          sstr >> sparsity;
          sstr >> time_limit_for_fix;
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

            std::cerr << std::endl << "TC #" << i << "..." << std::endl;
            
            std::vector<Point> points;
          

            //points =
              //generate_points_on_circle_2<Kernel>(NUM_POINTS, 3.);
              //generate_points_on_moment_curve<Kernel>(NUM_POINTS, ambient_dim, 0., 1.);
              //generate_points_on_plane<Kernel>(NUM_POINTS);
              //generate_points_on_sphere_3<Kernel>(NUM_POINTS, 3.0);
              //generate_points_on_sphere_d<Kernel>(NUM_POINTS, ambient_dim, 3.0);
              //generate_points_on_klein_bottle_3D<Kernel>(NUM_POINTS, 4., 3.);
              //generate_points_on_klein_bottle_4D<Kernel>(NUM_POINTS, 4., 3., NOISE);
              //generate_points_on_klein_bottle_variant_5D<Kernel>(NUM_POINTS, 4., 3.);

            /*if (input == "Klein_function")
              make_mesh_implicit(facet_approx, facet_sizing, cell_sizing, Klein_function(), input);
            else*/
            {
              load_points_from_file<Point>(
                input, std::back_inserter(points)/*, 600*/);
              make_tc(points, intrinsic_dim, sparsity, time_limit_for_fix);
            }

            std::cerr << "TC #" << i++ << " done." << std::endl;
            std::cerr << std::endl << "---------------------------------" << std::endl << std::endl;

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
    std::cerr << "Script file '" << BENCHMARK_SCRIPT_FILENAME << "' NOT found." << std::endl;
  }

  return 0;
}
