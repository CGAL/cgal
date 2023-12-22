//#define CHECK_MEMORY_LEAKS_ON_MSVC
#if defined(CHECK_MEMORY_LEAKS_ON_MSVC) && defined(_MSC_VER)
  #define _CRTDBG_MAP_ALLOC
  #include <stdlib.h>
  #include <crtdbg.h>
#endif

// Without TBB_USE_THREADING_TOOL Intel Inspector XE will report false positives in Intel TBB
// (https://www.intel.com/content/www/us/en/developer/articles/technical/compiler-settings-for-threading-error-analysis-in-intel-inspector-xe.html)
#ifdef _DEBUG
# define TBB_USE_THREADING_TOOL
#endif

#include <CGAL/Mesh_3/config.h>
#include <CGAL/Memory_sizer.h>

// ==========================================================================
// BENCHMARK GENERAL PARAMETERS
// ==========================================================================

// #define CGAL_MESH_3_BENCHMARK_EXPORT_TO_MAYA
#define CGAL_MESH_3_BENCHMARK_EXPORT_TO_MESH

// ==========================================================================
// MESH_3 GENERAL PARAMETERS
// ==========================================================================

// #define CGAL_MESH_3_POLYHEDRON_WITH_FEATURES
// #define CGAL_MESH_3_IMPLICIT_WITH_FEATURES

// #define CGAL_MESH_3_USE_OLD_SURFACE_RESTRICTED_DELAUNAY_UPDATE // WARNING: VERY SLOW
#define CGAL_MESH_3_INITIAL_POINTS_NO_RANDOM_SHOOTING

// #define CGAL_MESH_3_BENCHMARK_LLOYD
// #define CGAL_MESH_3_BENCHMARK_ODT
// #define CGAL_MESH_3_BENCHMARK_PERTURB
// #define CGAL_MESH_3_BENCHMARK_EXUDE
// #define CGAL_MESH_3_MANIFOLD

// #define CGAL_MESH_3_VERBOSE
// #define CGAL_MESH_3_PERTURBER_VERBOSE
// #define CGAL_MESH_3_PERTURBER_HIGH_VERBOSITY
// #define CGAL_MESH_3_EXUDER_VERBOSE
// #define CGAL_MESH_3_EXUDER_HIGH_VERBOSITY
// #define CGAL_MESH_3_VERY_VERBOSE
// #define CGAL_MESHES_DEBUG_REFINEMENT_POINTS
// #define CGAL_MESH_3_OPTIMIZER_VERBOSE

#define CGAL_MESH_3_PROFILING
// #define CHECK_AND_DISPLAY_THE_NUMBER_OF_BAD_ELEMENTS_IN_THE_END

// ==========================================================================
// INPUTS
// ==========================================================================

const char* const BENCHMARK_INPUTS_FILENAME = "benchmark_inputs.txt";

const bool USE_RELATIVE_CRITERIA_VALUES = true; // relative to the bbox's diagonal
const double DEFAULT_FACE_SIZE          = 100;   // can be relative
const double DEFAULT_FACE_APPROX        = 200;  // can be relative
const double DEFAULT_FACE_ANGLE         = 25;   // cannot be relative
const double DEFAULT_CELL_SIZE          = 100;   // can be relative
const double DEFAULT_CELL_RADIUS_RATIO  = 3;    // cannot be relative

// ==========================================================================
// CONCURRENCY
// ==========================================================================

#ifdef CGAL_CONCURRENT_MESH_3

# include <tbb/task_arena.h>

# ifndef CGAL_LINKED_WITH_TBB
#  pragma message(" : Warning: CGAL_LINKED_WITH_TBB not defined: EVERYTHING WILL BE SEQUENTIAL.")
# endif

// # define BENCHMARK_WITH_1_TO_MAX_THREADS

// # define CGAL_MESH_3_USE_LAZY_SORTED_REFINEMENT_QUEUE // default behavior
// # define CGAL_MESH_3_USE_LAZY_UNSORTED_REFINEMENT_QUEUE
// # define CGAL_MESH_3_IF_UNSORTED_QUEUE_JUST_SORT_AFTER_SCAN // recommended
// # define CGAL_PARALLEL_MESH_3_DO_NOT_ADD_OUTSIDE_POINTS_ON_A_FAR_SPHERE // not recommended

  // ==========================================================================
  // Verbose
  // ==========================================================================

# define CGAL_CONCURRENT_MESH_3_VERBOSE
# define CGAL_CONCURRENT_MESH_3_VERY_VERBOSE

  // ==========================================================================
  // Concurrency config
  // ==========================================================================

  const char* const CONCURRENT_MESHER_CONFIG_FILENAME = "concurrent_mesher_config.cfg";

  // =====================
  // Worksharing strategy
  // =====================

// #define CGAL_MESH_3_LOAD_BASED_WORKSHARING // Not recommended

  // ==========================================================================
  // Profiling
  // ==========================================================================

  // For profiling, etc.
# define CGAL_CONCURRENT_MESH_3_PROFILING
// # define CGAL_DEBUG_FORCE_SEQUENTIAL_MESH_REFINEMENT

#include <thread>

// ==========================================================================
// SEQUENTIAL
// ==========================================================================

#else // !CGAL_CONCURRENT_MESH_3

// # define CGAL_MESH_3_USE_LAZY_SORTED_REFINEMENT_QUEUE // default behavior
// # define CGAL_MESH_3_USE_LAZY_UNSORTED_REFINEMENT_QUEUE
// # define CGAL_MESH_3_IF_UNSORTED_QUEUE_JUST_SORT_AFTER_SCAN // recommended

#endif // CGAL_CONCURRENT_MESH_3