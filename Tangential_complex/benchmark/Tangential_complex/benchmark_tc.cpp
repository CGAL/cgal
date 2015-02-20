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

#include <boost/algorithm/string/replace.hpp>
#include <boost/algorithm/string/trim_all.hpp>

#include <cstdlib>
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
    subelements.push_back("Sparsity");
    subelements.push_back("Num_points_in_input");
    subelements.push_back("Num_points");
    subelements.push_back("Initial_num_inconsistent_local_tr");
    subelements.push_back("Best_num_inconsistent_local_tr");
    subelements.push_back("Final_num_inconsistent_local_tr");
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

void make_tc(std::vector<Point> &points, int intrinsic_dim,
             double sparsity = 0., double time_limit_for_fix = 0.,
             const char *input_name = "tc")
{
  Kernel k;
  Wall_clock_timer t;

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
    
#ifdef CGAL_TC_PROFILING
  Wall_clock_timer t_gen;
#endif
  
#ifdef CGAL_TC_PROFILING
  std::cerr << "Point set generated in " << t_gen.elapsed()
            << " seconds." << std::endl;
#endif
  
  CGAL_TC_SET_PERFORMANCE_DATA("Num_points_in_input", points.size());
  
#ifdef USE_ANOTHER_POINT_SET_FOR_TANGENT_SPACE_ESTIM
  std::vector<Point> points_not_sparse = points;
#endif

  if (sparsity != 0.)
  {
    std::size_t num_points_before = points.size();
    points = sparsify_point_set(k, points, sparsity*sparsity);
    std::cerr << "Number of points before/after sparsification: "
      << num_points_before << " / " << points.size() << std::endl;
  }
  
  CGAL_TC_SET_PERFORMANCE_DATA("Sparsity", sparsity);
  CGAL_TC_SET_PERFORMANCE_DATA("Num_points", points.size());
  
#ifdef USE_ANOTHER_POINT_SET_FOR_TANGENT_SPACE_ESTIM
  TC tc(points.begin(), points.end(), sparsity, intrinsic_dim,
    points_not_sparse.begin(), points_not_sparse.end(), k);
#else
  TC tc(points.begin(), points.end(), sparsity, intrinsic_dim, k);
#endif

  double init_time = t.elapsed(); t.reset();

  tc.compute_tangential_complex();
  double computation_time = t.elapsed(); t.reset();
    
  //tc.check_if_all_simplices_are_in_the_ambient_delaunay();

  double export_before_time = -1.;
  if (intrinsic_dim <= 3)
  {
    t.reset();
    std::stringstream output_filename;
    output_filename << "output/" << input_name_stripped << "_" << intrinsic_dim
      << "_in_R" << ambient_dim << "_BEFORE_FIX.off";
    std::ofstream off_stream(output_filename.str().c_str());
    tc.export_to_off(off_stream, true);
    export_before_time = t.elapsed(); t.reset();
  }


  t.reset();
  unsigned int num_fix_steps;
  std::size_t initial_num_inconsistent_local_tr;
  std::size_t best_num_inconsistent_local_tr;
  std::size_t final_num_inconsistent_local_tr;
  CGAL::Fix_inconsistencies_status fix_ret = tc.fix_inconsistencies(
    num_fix_steps, initial_num_inconsistent_local_tr,
    best_num_inconsistent_local_tr, final_num_inconsistent_local_tr,
    time_limit_for_fix);
  double fix_time = t.elapsed(); t.reset();

  CGAL_TC_SET_PERFORMANCE_DATA("Initial_num_inconsistent_local_tr", 
                               initial_num_inconsistent_local_tr);
  CGAL_TC_SET_PERFORMANCE_DATA("Best_num_inconsistent_local_tr", 
                               best_num_inconsistent_local_tr);
  CGAL_TC_SET_PERFORMANCE_DATA("Final_num_inconsistent_local_tr", 
                               final_num_inconsistent_local_tr);
  
  tc.check_and_solve_inconsistencies_by_adding_higher_dim_simplices();
  TC::Simplicial_complex complex;
  int max_dim = tc.export_TC(complex);
  tc.check_if_all_simplices_are_in_the_ambient_delaunay(&complex);
  complex.collapse(max_dim);
  complex.display_stats();

  std::ofstream off_stream("output/test.off"); // CJTODO TEMP TEST
  tc.export_to_off(complex, off_stream);
  
  double export_after_time = -1.;
  if (intrinsic_dim <= 3)
  {
    t.reset();
    std::stringstream output_filename;
    output_filename << "output/" << input_name_stripped << "_" << intrinsic_dim
      << "_in_R" << ambient_dim << "_AFTER_FIX.off";
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
          boost::replace_all(line, "\t", " ");
          boost::trim_all(line);
          std::cerr << std::endl << std::endl;
          std::cerr << "*****************************************" << std::endl;
          std::cerr << "******* " << line << std::endl;
          std::cerr << "*****************************************" << std::endl;
          std::stringstream sstr(line);

          std::string input;
          std::string param1;
          std::string param2;
          std::string param3;
          std::size_t num_points;
          int ambient_dim;
          int intrinsic_dim;
          double sparsity;
          double time_limit_for_fix;
          int num_iteration;
          sstr >> input;
          sstr >> param1;
          sstr >> param2;
          sstr >> param3;
          sstr >> num_points;
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

            if (input == "generate_moment_curve")
            {
              points = generate_points_on_moment_curve<Kernel>(
                num_points, ambient_dim, 
                std::atof(param1.c_str()), std::atof(param2.c_str()));
            }
            else if (input == "generate_plane")
            {
              points = generate_points_on_plane<Kernel>(num_points);
            }
            else if (input == "generate_sphere_d")
            {
              points = generate_points_on_sphere_d<Kernel>(
                num_points, ambient_dim, 
                std::atof(param1.c_str()));
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
              load_points_from_file<Point>(
                input, std::back_inserter(points)/*, 600*/);
            }

            if (!points.empty())
            {
              make_tc(points, intrinsic_dim, sparsity, 
                      time_limit_for_fix, input.c_str());

              std::cerr << "TC #" << i++ << " done." << std::endl;
              std::cerr << std::endl << "---------------------------------" 
                        << std::endl << std::endl;
            }
            else
            {
              std::cerr << "TC #" << i++ << ": no points loaded." << std::endl;
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
    std::cerr << "Script file '" << BENCHMARK_SCRIPT_FILENAME << "' NOT found." << std::endl;
  }

  return 0;
}
