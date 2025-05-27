#include "benchmark_config.h"
#include "benchmark_xml.h"
#include "mesh_quality.h"

std::string XML_perf_data::default_filename;

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/tetrahedral_remeshing.h>
#include <CGAL/Tetrahedral_remeshing/Remeshing_triangulation_3.h>

#include <CGAL/Real_timer.h>

#include <CGAL/boost/graph/IO/polygon_mesh_io.h>
#include <CGAL/Polygon_mesh_processing/self_intersections.h>


#ifdef CGAL_LINKED_WITH_TBB
# include <tbb/task_arena.h>
# define TBB_PREVIEW_GLOBAL_CONTROL 1
# include <tbb/global_control.h>
#endif


#include <filesystem>
#include <fstream>
#include <iostream>
#include <signal.h>
#include <string>
#include <sstream>

namespace PMP = CGAL::Polygon_mesh_processing;

// basic types from kernel
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::FT FT;
typedef K::Point_3 Point;
typedef K::Sphere_3 Sphere;

typedef CGAL::Tetrahedral_remeshing::Remeshing_triangulation_3<K> Remeshing_triangulation;
// To avoid verbose function and named parameters call
using namespace CGAL::parameters;

#include "implicit_functions.h"

std::string get_output_filename(const std::string& input_name)
{
  std::string filename = std::string(input_name);
  filename = filename.substr(filename.find_last_of("/") + 1, filename.length() - 1);
  filename = filename.substr(0, filename.find_last_of("."));
  return filename;
}

std::string get_technique()
{
  // This function should be implemented when we implement our parallel framework
  return "Sequential remeshing";
}

void display_info(int num_threads)
{
  // This function should be implemented when we implement our parallel framework
  std::string tech= get_technique();
  std::cout << tech << std::endl;
}

void xml_perf_set_technique()
{
  CGAL_TETRAHEDRAL_REMESHING_SET_PERFORMANCE_DATA("Technique", get_technique());
}

enum Exit_code
{
  // Success
  VALID_OUTPUT = 0,

  // Failure
  INPUT_IS_INVALID = 1,
  OUTPUT_IS_INVALID = 2
};

#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Polyhedral_mesh_domain_with_features_3.h>
Exit_code remesh(const std::string& input_filename,
                                          double edge_length,int num_iterations)
{

    Remeshing_triangulation tr;
    std::ifstream is(input_filename, std::ios_base::in);
    if(!CGAL::IO::read_MEDIT(is,tr)){
      std::cerr << "Error: Could not read '" << input_filename  << std::endl;
      return INPUT_IS_INVALID;
    }
  
  generate_quality_metrics(tr);
  XML_perf_data::commit();
  
  std::cout << "remeshing:" << input_filename  << std::endl;

  CGAL::Real_timer t;
  t.start();
    CGAL::tetrahedral_isotropic_remeshing(tr,edge_length,number_of_iterations(num_iterations));
  t.stop();
  std::cout << "Remeshing done." << std::endl;

  CGAL_TETRAHEDRAL_REMESHING_SET_PERFORMANCE_DATA("V", tr.number_of_vertices());
  CGAL_TETRAHEDRAL_REMESHING_SET_PERFORMANCE_DATA("F", tr.number_of_facets());
  CGAL_TETRAHEDRAL_REMESHING_SET_PERFORMANCE_DATA("C",tr.number_of_cells());
  CGAL_TETRAHEDRAL_REMESHING_SET_PERFORMANCE_DATA("Mem", CGAL::Memory_sizer().virtual_size() >> 20);
  CGAL_TETRAHEDRAL_REMESHING_SET_PERFORMANCE_DATA("Total_time", t.time());

#ifdef CGAL_TETRAHEDRAL_REMESHING_BENCHMARK_EXPORT_TO_MESH
  const std::string& output_filename=get_output_filename(input);
  std::cout << "Exporting to " << output_filename + ".mesh (Medit)... ";
  std::ofstream out_medit(output_filename + ".mesh");
  CGAL::IO::write_MEDIT(out_medit, tr);
  std::cout << "done." << std::endl;
#endif

generate_quality_metrics(tr);

XML_perf_data::commit();
return VALID_OUTPUT;
}

Exit_code run_tetrahedral_remeshing(const std::string& input,
                            double target_edge_length,int num_iterations)
                            {
    std::string domain = input;
    size_t slash_index = domain.find_last_of('/');
    if(slash_index == std::string::npos)
      slash_index = domain.find_last_of('\\');
    if(slash_index == std::string::npos)
      slash_index = 0;
    else
      ++slash_index;

    domain = domain.substr(slash_index, domain.find_last_of('.') - slash_index);

    CGAL_TETRAHEDRAL_REMESHING_SET_PERFORMANCE_DATA("Domain", domain);
    CGAL_TETRAHEDRAL_REMESHING_SET_PERFORMANCE_DATA("target_edge_length",target_edge_length);
    CGAL_TETRAHEDRAL_REMESHING_SET_PERFORMANCE_DATA("num_iterations",num_iterations);
    xml_perf_set_technique();
    Exit_code res;
    res=remesh(input,target_edge_length,num_iterations);
    return res;                  
  }

// reads filename & input parameters either from argv or from a script file
int main(int argc, char** argv)
{

  std::cout.precision(17);
  std::cerr.precision(17);
  // for the default xml filename, check XML_perf_data::build_filename()
  if(argc > 1)
    XML_perf_data::default_filename = argv[1];

#if defined(CHECK_MEMORY_LEAKS_ON_MSVC) && defined(_MSC_VER)
  _CrtSetDbgFlag ( _CRTDtbbBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF );
#endif

#ifdef CGAL_CONCURRENT_TETRAHEDRAL_REMESHING
  Concurrent_tetrahedral_remesher_config::load_config_file(CONCURRENT_TETRAHEDRAL_REMESHING_CONFIG_FILENAME, true);

  int max_nb_threads = Concurrent_tetrahedral_remeshing_config::get().num_threads;
  if(max_nb_threads == -1) // if not set in the config file, take the max available
    max_nb_threads = tbb::this_task_arena::max_concurrency();
#endif

  Exit_code script_res = VALID_OUTPUT;

#ifdef CGAL_CONCURRENT_TETRAHEDRAL_REMESHING
 #ifdef BENCHMARK_WITH_1_TO_MAX_THREADS
  for(int num_threads=1; num_threads<=max_nb_threads; ++num_threads)
 #else
  int num_threads = max_nb_threads;
 #endif // BENCHMARK_WITH_1_TO_MAX_THREADS
#endif // CGAL_CONCURRENT_TETRAHEDRAL_REMESHING
  {
#ifdef CGAL_CONCURRENT_TETRAHEDRAL_REMESHING
    std::cout << "-- Parallel Tetrahedral Remeshing --" << std::endl;
    tbb::global_control control(tbb::global_control::max_allowed_parallelism, num_threads);
    display_info(num_threads);

    CGAL_TETRAHEDRAL_REMESHING_SET_PERFORMANCE_DATA("Num_threads", num_threads);
    CGAL_TETRAHEDRAL_REMESHING_SET_PERFORMANCE_DATA("Lockgrid_size", Concurrent_tetrahedral_remesher_config::get().locking_grid_num_cells_per_axis);
    CGAL_TETRAHEDRAL_REMESHING_SET_PERFORMANCE_DATA("Lock_radius", Concurrent_tetrahedral_remesher_config::get().first_grid_lock_radius);
    CGAL_TETRAHEDRAL_REMESHING_SET_PERFORMANCE_DATA("Num_work_items_per_batch", Concurrent_tetrahedral_remesher_config::get().num_work_items_per_batch);
#else
    std::cout << "-- Sequential Tetrahedral Remeshing --" << std::endl;
    
    CGAL_TETRAHEDRAL_REMESHING_SET_PERFORMANCE_DATA("Num_threads", "N/A");
    CGAL_TETRAHEDRAL_REMESHING_SET_PERFORMANCE_DATA("Lockgrid_size", "N/A");
    CGAL_TETRAHEDRAL_REMESHING_SET_PERFORMANCE_DATA("Lock_radius", "N/A");
    CGAL_TETRAHEDRAL_REMESHING_SET_PERFORMANCE_DATA("Statgrid_size", "N/A");
    CGAL_TETRAHEDRAL_REMESHING_SET_PERFORMANCE_DATA("Num_work_items_per_batch", "N/A");
#endif

    // Script file format: each line gives
    //    - Filename (polyhedron and image) or "XXX_function" (implicit)
    //    - Facet and cell criteria
    //    - Number of iterations with these parameters
    std::ifstream script_file;
    script_file.open(BENCHMARK_INPUTS_FILENAME);
    if(script_file.is_open())
    {
      std::cout << "Found inputs file '" << argv[1]<< "'" << std::endl;

      std::string line;
      while(std::getline(script_file, line))
      {
        if(line.empty() || line[0] == '#') // lines starting with '#' are ignored
          continue;

        std::cout << std::endl << std::endl;
        std::cout << "*****************************************" << std::endl;
        std::cout << "******* " << line << std::endl;
        std::cout << "*****************************************" << std::endl;

        std::stringstream sstr(line);

        std::string input;
        double target_edge_length;
        int num_iterations;
        if(!(sstr >> input
                  >> target_edge_length
                  >> num_iterations
                  ))
        {
          std::cerr << "Error: failed to read input" << std::endl;
          return INPUT_IS_INVALID;
        }

        script_res = run_tetrahedral_remeshing(input, target_edge_length,
                                               num_iterations);
      }
    }
    else // no script
    {
      std::cout << "Inputs file '" << BENCHMARK_INPUTS_FILENAME << "' NOT found." << std::endl;

      // If the script is empty, use the command line arguments:
      // [this_program]
      // - xml_filename
      // - filename
      // - target_edge_length
      // - num_iterations

      // Custom: Use command line arguments for Tetrahedral_remeshing
      std::string input = (argc > 2) ? argv[2] : "";
      double target_edge_length = (argc > 3) ? std::stod(argv[3]) : 0.1;
      int num_iterations = (argc > 4) ? std::stoi(argv[4]) : 1;

      if (input.empty()) {
        std::cerr << "Error: No input mesh file provided." << std::endl;
        return INPUT_IS_INVALID;
      }

      script_res = run_tetrahedral_remeshing(input, target_edge_length, num_iterations);
    }
    
  }

  return script_res;
}
