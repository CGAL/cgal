
#include <string>

#include <boost/program_options.hpp>
namespace po = boost::program_options;

const char * const BENCHMARK_CONFIG_FILENAME = 
  "D:/INRIA/CGAL/workingcopy/Mesh_3/benchmark/Mesh_3/concurrency_config.cfg";


// ==========================================================================
// MESH_3 GENERAL PARAMETERS
// ==========================================================================

//#define CGAL_MESH_3_VERBOSE
//#define CGAL_MESH_3_VERY_VERBOSE
#define CGAL_MESH_3_LAZY_REFINEMENT_QUEUE
#define CGAL_MESH_3_INITIAL_POINTS_NO_RANDOM_SHOOTING

// ==========================================================================
// CONCURRENCY
// ==========================================================================

#ifdef CONCURRENT_MESH_3

# include <tbb/tbb.h>

  // ==========================================================================
  // Concurrency activation
  // ==========================================================================

# define CGAL_MESH_3_CONCURRENT_SCAN_TRIANGULATION
# define CGAL_MESH_3_CONCURRENT_REFINEMENT
  // In case some code uses CGAL_PROFILE, it needs to be concurrent
# define CGAL_CONCURRENT_PROFILE
# define CGAL_CONCURRENT_MESH_3_VERBOSE
//#define CGAL_CONCURRENT_MESH_3_VERY_VERBOSE

  // ==========================================================================
  // Locking strategy
  // ==========================================================================

# ifdef CGAL_MESH_3_CONCURRENT_REFINEMENT

    const char * const CONFIG_FILENAME = 
      "D:/INRIA/CGAL/workingcopy/Mesh_3/demo/Mesh_3/concurrent_mesher_config.cfg";
    
//#   define CGAL_MESH_3_LOCKING_STRATEGY_CELL_LOCK
#   define CGAL_MESH_3_LOCKING_STRATEGY_SIMPLE_GRID_LOCKING

//#   define CGAL_MESH_3_CONCURRENT_REFINEMENT_LOCK_ADJ_CELLS // USELESS, FOR TESTS ONLY
//#   define CGAL_MESH_3_DO_NOT_LOCK_INFINITE_VERTEX // DOES NOT WORK YET
//#   define CGAL_MESH_3_ACTIVATE_GRID_INDEX_CACHE_IN_VERTEX // DOES NOT WORK YET

//#   define CGAL_MESH_3_WORKSHARING_USES_TASKS
//#   define CGAL_MESH_3_WORKSHARING_USES_PARALLEL_FOR
#   define CGAL_MESH_3_WORKSHARING_USES_PARALLEL_DO
    
# define CGAL_MESH_3_ADD_OUTSIDE_POINTS_ON_A_FAR_SPHERE

#   ifdef CGAL_MESH_3_LOCKING_STRATEGY_CELL_LOCK
#     include <tbb/recursive_mutex.h>
      typedef tbb::recursive_mutex Cell_mutex_type; // CJTODO try others
#   endif

# endif

  // ==========================================================================
  // CJTODO TEMP
  // ==========================================================================
# include <tbb/tbb.h>
  typedef tbb::recursive_mutex Global_mutex_type;
  extern Global_mutex_type g_global_mutex; // CJTODO: temporary

  // ==========================================================================
  // Profiling
  // ==========================================================================

  // For abortion profiling, etc.
# define CGAL_CONCURRENT_MESH_3_PROFILING
  
  // ==========================================================================
  // TBB
  // ==========================================================================

  // Use TBB malloc proxy (for all new/delete/malloc/free calls)
# include <tbb/tbbmalloc_proxy.h>
  

#else // !CONCURRENT_MESH_3

# define CGAL_MESH_3_USE_LAZY_UNSORTED_REFINEMENT_QUEUE
# define CGAL_MESH_3_IF_UNSORTED_QUEUE_JUST_SORT_AFTER_SCAN

#endif // CONCURRENT_MESH_3
  

#define MESH_3_PROFILING
#define CHECK_AND_DISPLAY_THE_NUMBER_OF_BAD_ELEMENTS_IN_THE_END
  
// ==========================================================================
// ==========================================================================
  
const char * const DEFAULT_INPUT_FILE_NAME = "D:/INRIA/CGAL/workingcopy/Mesh_3/examples/Mesh_3/data/elephant.off";
//const char *DEFAULT_INPUT_FILE_NAME = "D:/INRIA/CGAL/workingcopy/Mesh_3/examples/Mesh_3/data/fandisk.off";
  
#ifdef CONCURRENT_MESH_3
  // CJTODO TEMP TEST
#ifdef CGAL_MESH_3_DO_NOT_LOCK_INFINITE_VERTEX
  bool g_is_set_cell_active = true;
#endif

  Global_mutex_type g_global_mutex; // CJTODO: temporary
  
  // CJTODO TEMP: not thread-safe => move it to Mesher_3
  
# ifdef CGAL_MESH_3_LOCKING_STRATEGY_CELL_LOCK
#   include <utility>
#   include <vector>
#   include <tbb/enumerable_thread_specific.h>
    tbb::enumerable_thread_specific<std::vector<std::pair<void*, unsigned int> > > g_tls_locked_cells;
# endif

#endif


// ==========================================================================
// ==========================================================================



#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>

#include <CGAL/Polyhedral_mesh_domain_3.h>
#include <CGAL/Implicit_mesh_domain_3.h>
#include <CGAL/Mesh_domain_with_polyline_features_3.h>

#include <CGAL/make_mesh_3.h>
#include <CGAL/refine_mesh_3.h>
    
// basic types from kernel
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::FT FT;
typedef Kernel::Point_3 Point;
typedef Kernel::Sphere_3 Sphere;
// IO
#include <CGAL/IO/Polyhedron_iostream.h>

// To avoid verbose function and named parameters call
using namespace CGAL::parameters;

struct Mesh_parameters
{
  double facet_angle;
  double facet_sizing;
  double facet_approx;
  
  double tet_shape;
  double tet_sizing;
  
  std::string log() const
  {
    std::stringstream sstr;
    sstr
      << " * facet min angle: " << facet_angle << std::endl
      << " * facet max size: " << facet_sizing << std::endl
      << " * facet approx error: " << facet_approx << std::endl
      << " * tet shape (radius-edge): " << tet_shape << std::endl
      << " * tet max size: " << tet_sizing << std::endl;

    return sstr.str();
  }
};


struct Klein_function
{
  typedef ::FT           FT;
  typedef ::Point        Point;
  
  FT operator()(const Point& query) const
  { 
	  const FT x = query.x();
	  const FT y = query.y();
	  const FT z = query.z();

    return   (x*x+y*y+z*z+2*y-1)
           * ( (x*x+y*y+z*z-2*y-1) *(x*x+y*y+z*z-2*y-1)-8*z*z)
           + 16*x*z* (x*x+y*y+z*z-2*y-1);
  }
};


bool make_mesh_polyhedron(const std::string &input_filename, 
                 double facet_sizing, 
                 double cell_sizing)
{
  // Domain
  typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
  typedef CGAL::Polyhedron_3<K> Polyhedron;
  typedef CGAL::Polyhedral_mesh_domain_3<Polyhedron, K> Mesh_domain;

  // Triangulation
  typedef CGAL::Mesh_triangulation_3<Mesh_domain>::type Tr;

  typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr> C3t3;

  // Criteria
  typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;

  // Create input polyhedron
  Polyhedron polyhedron;
  std::ifstream input(input_filename);
  if (!input.is_open())
  {
    std::cerr << "Could not open '" << input_filename << "'" << std::endl;
    return false;
  }
  input >> polyhedron;
   
  // Create domain
  Mesh_domain domain(polyhedron);

  Mesh_parameters params;
  params.facet_angle = 25;
  params.facet_sizing = facet_sizing;
  params.facet_approx = 0.0068;
  params.tet_sizing = cell_sizing;
  params.tet_shape = 3;

  std::cerr 
    << "File: " << input_filename << std::endl
    << "Parameters: " << std::endl 
    << params.log() << std::endl;

  // Mesh criteria (no cell_size set)
  Mesh_criteria criteria(
    facet_angle=params.facet_angle,
    facet_size=params.facet_sizing,
    facet_distance=params.facet_approx,
    cell_size=params.tet_sizing,
    cell_radius_edge_ratio=params.tet_shape
  );

  // Mesh generation
  C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria, no_perturb(), no_exude());

  std::cerr
    << "Vertices: " << c3t3.triangulation().number_of_vertices() << std::endl
    << "Facets  : " << c3t3.number_of_facets_in_complex() << std::endl
    << "Tets    : " << c3t3.number_of_cells_in_complex() << std::endl;

  return true;
}

bool make_mesh_implicit(double facet_sizing, double cell_sizing)
{
  // Domain
  typedef CGAL::Implicit_mesh_domain_3<const Klein_function, Kernel> Mesh_domain;
  // Triangulation
  typedef CGAL::Mesh_triangulation_3<Mesh_domain>::type Tr;
  typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr> C3t3;
  // Criteria
  typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;
     
  // Create domain
  Klein_function f;
	Sphere bounding_sphere(CGAL::ORIGIN, 7.0 * 7.0);
  Mesh_domain domain(f, bounding_sphere);

  Mesh_parameters params;
  params.facet_angle = 25;
  params.facet_sizing = facet_sizing;
  params.facet_approx = 0.0068;
  params.tet_sizing = cell_sizing;
  params.tet_shape = 3;
  
  std::cerr 
    << "Klein function" << std::endl
    << "Parameters: " << std::endl 
    << params.log() << std::endl;

  // Mesh criteria (no cell_size set)
  Mesh_criteria criteria(
    facet_angle=params.facet_angle,
    facet_size=params.facet_sizing,
    facet_distance=params.facet_approx,
    cell_size=params.tet_sizing,
    cell_radius_edge_ratio=params.tet_shape
  );

  // Mesh generation
  C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria, no_perturb(), no_exude());

  std::cerr
    << "Vertices: " << c3t3.triangulation().number_of_vertices() << std::endl
    << "Facets  : " << c3t3.number_of_facets_in_complex() << std::endl
    << "Tets    : " << c3t3.number_of_cells_in_complex() << std::endl;

  return true;
}

int main()
{
  // Program options
  po::variables_map vm;
  try
  {
    // Declare the supported options.
    po::options_description desc("Allowed options");
    desc.add_options()
      ("filename", po::value<std::string>()->default_value(DEFAULT_INPUT_FILE_NAME), "")
      ("facet_sizing", po::value<double>()->default_value(0.005), "")
      ("cell_sizing", po::value<double>()->default_value(0.005), "")
      ("numthreads", po::value<int>()->default_value(-1), "");

    po::store(po::parse_config_file<char>(BENCHMARK_CONFIG_FILENAME, desc), vm);
    po::notify(vm); 
  }
  catch (std::exception &e)
  {
    std::cerr << "Config file error: " << e.what() << std::endl;
    return false;
  }
  int num_threads = vm["numthreads"].as<int>();
  double facet_sizing = vm["facet_sizing"].as<double>();
  double cell_sizing = vm["cell_sizing"].as<double>();
  std::string filename = vm["filename"].as<std::string>();

#ifdef CONCURRENT_MESH_3
  tbb::task_scheduler_init init(tbb::task_scheduler_init::deferred);
  if (num_threads > 0)
    init.initialize(num_threads);
#endif

  for(int i = 1 ; ; ++i)
  {
    std::cerr << "Refinement #" << i << "..." << std::endl;
    
#ifdef CONCURRENT_MESH_3
    
    std::cerr << "CONCURRENT MESH_3" << std::endl;

# if defined(CGAL_MESH_3_WORKSHARING_USES_TASKS)
    std::cerr << "Using TBB task-scheduler" << std::endl;
# elif defined(CGAL_MESH_3_WORKSHARING_USES_PARALLEL_FOR)
    std::cerr << "Using tbb::parallel_for" << std::endl;
# elif defined(CGAL_MESH_3_WORKSHARING_USES_PARALLEL_DO)
    std::cerr << "Using tbb::parallel_do" << std::endl;
# else
    std::cerr << "Using unknown technique" << std::endl;
# endif
    if (num_threads != -1)
      std::cerr << "Num threads = " << num_threads << std::endl;
    else
      std::cerr << "Num threads = AUTO" << std::endl;

#else
    std::cerr << "SEQUENTIAL MESH_3";
# if defined(CGAL_MESH_3_USE_LAZY_UNSORTED_REFINEMENT_QUEUE)
    std::cerr << "(unsorted refinement queue, ";
# else
    std::cerr << "(sorted refinement queue, ";
# endif
# if defined(CGAL_MESH_3_INITIAL_POINTS_NO_RANDOM_SHOOTING)
    std::cerr << "NO random shooting)" << std::endl;
# else
    std::cerr << "WITH random shooting)" << std::endl;
# endif
    
#endif

    //make_mesh_polyhedron(filename, facet_sizing, cell_sizing);
    make_mesh_implicit(facet_sizing, cell_sizing);
    std::cerr << "Refinement #" << i << " done." << std::endl;
    std::cerr << std::endl << "---------------------------------" << std::endl << std::endl;
  }

  return 0;
}
