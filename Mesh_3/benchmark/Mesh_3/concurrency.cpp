
// ==========================================================================
// CONCURRENCY
// ==========================================================================

#ifdef CONCURRENT_MESH_3

  // ==========================================================================
  // Concurrency activation
  // ==========================================================================

# define CGAL_MESH_3_LAZY_REFINEMENT_QUEUE
# define CGAL_MESH_3_CONCURRENT_SCAN_TRIANGULATION
# define CGAL_MESH_3_CONCURRENT_REFINEMENT
  // In case some code uses CGAL_PROFILE, it needs to be concurrent
# define CGAL_CONCURRENT_PROFILE
//# define CGAL_CONCURRENT_MESH_3_VERBOSE

  // ==========================================================================
  // Locking strategy
  // ==========================================================================

# ifdef CGAL_MESH_3_CONCURRENT_REFINEMENT
//#   define CGAL_MESH_3_LOCKING_STRATEGY_CELL_LOCK
#   define CGAL_MESH_3_LOCKING_STRATEGY_SIMPLE_GRID_LOCKING
//#   define CGAL_MESH_3_CONCURRENT_REFINEMENT_LOCK_ADJ_CELLS

    const int MESH_3_LOCKING_GRID_NUM_CELLS_PER_AXIS = 30;
    const int MESH_3_FIRST_GRID_LOCK_RADIUS = 2;
    const int MESH_3_REFINEMENT_GRAINSIZE = 10;

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
  // Concurrency Parameters
  // ==========================================================================

  const size_t ELEMENT_BATCH_SIZE = 10000;

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

#endif // CONCURRENT_MESH_3
  
#define MESH_3_PROFILING

  
// ==========================================================================
// ==========================================================================
  
// CJTODO TEMP
bool g_temp = false;

#ifdef CONCURRENT_MESH_3
  #include <CGAL/Mesh_3/Locking_data_structures.h> // CJODO TEMP?
  #include <CGAL/BBox_3.h>

  // CJTODO TEMP TEST
#ifdef CGAL_MESH_3_DO_NOT_LOCK_INFINITE_VERTEX
  bool g_is_set_cell_active = true;
#endif

  Global_mutex_type g_global_mutex; // CJTODO: temporary
  
  // CJTODO TEMP: not thread-safe => move it to Mesher_3
  // Elephant.off => BBox (x,y,z): [ -0.358688, 0.356308 ], [ -0.498433, 0.49535 ], [ -0.298931, 0.298456 ]
  CGAL::Bbox_3 g_bbox(-0.35, 0.35, -0.5, 0.5, -0.3, 0.3);
# ifdef CGAL_MESH_3_LOCKING_STRATEGY_SIMPLE_GRID_LOCKING
  CGAL::Mesh_3::Refinement_grid_type g_lock_grid(g_bbox, MESH_3_LOCKING_GRID_NUM_CELLS_PER_AXIS);

# elif defined(CGAL_MESH_3_LOCKING_STRATEGY_CELL_LOCK)
# include <utility>
# include <vector>
# include <tbb/enumerable_thread_specific.h>
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
#include <CGAL/make_mesh_3.h>
#include <CGAL/refine_mesh_3.h>

// IO
#include <CGAL/IO/Polyhedron_iostream.h>

// Domain
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Polyhedron_3<K> Polyhedron;
typedef CGAL::Polyhedral_mesh_domain_3<Polyhedron, K> Mesh_domain;

// Triangulation
typedef CGAL::Mesh_triangulation_3<Mesh_domain>::type Tr;

typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr> C3t3;

// Criteria
typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;

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

bool refine_mesh(const std::string &input_filename)
{
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
  params.facet_sizing = 0.002;
  params.facet_approx = 0.0068;
  /*params.tet_shape = 3;
  params.tet_sizing = 1.;*/
  
  std::cerr 
    << "Parameters: " << std::endl 
    << params.log() << std::endl;

  // Mesh criteria (no cell_size set)
  Mesh_criteria criteria(
    facet_angle=params.facet_angle,
    facet_size=params.facet_sizing,
    facet_distance=params.facet_approx/*,
    cell_size=params.tet_sizing,
    cell_radius_edge_ratio=params.tet_shape*/
  );

  // Mesh generation
  C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria, no_perturb(), no_exude());

  // Output
  /*std::ofstream medit_file("out_1.mesh");
  c3t3.output_to_medit(medit_file);
  medit_file.close();*/

  /*
  // Set tetrahedron size (keep cell_radius_edge_ratio), ignore facets
  Mesh_criteria new_criteria(cell_radius_edge_ratio=3, cell_size=0.03);

  // Mesh refinement
  CGAL::refine_mesh_3(c3t3, domain, new_criteria);

  // Output
  medit_file.open("out_2.mesh");
  c3t3.output_to_medit(medit_file);*/

  return true;
}

int main()
{
  for(int i = 1 ; ; ++i)
  {
    std::cerr << "Refinement #" << i << "..." << std::endl;
    refine_mesh("D:/INRIA/CGAL/workingcopy/Mesh_3/examples/Mesh_3/data/elephant.off");
    std::cerr << "Refinement #" << i << " done." << std::endl;
    std::cerr << std::endl << "---------------------------------" << std::endl << std::endl;
  }

  return 0;
}
