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
#define CGAL_TETRAHEDRAL_REMESHING_BENCHMARK_EXPORT_TO_MESH

// ==========================================================================
// TETRAHEDRAL REMESHING GENERAL PARAMETERS
// ==========================================================================
#define CGAL_TETRAHEDRAL_REMESHING_VERBOSE
#define CGAL_TETRAHEDRAL_REMESHING_PROFILING

// ==========================================================================
// INPUTS
// ==========================================================================

const char* const BENCHMARK_INPUTS_FILENAME = "benchmark_inputs.txt";

// ==========================================================================
// CONCURRENCY
// ==========================================================================

#ifdef CGAL_CONCURRENT_TETRAHEDRAL_REMESHING

# include <tbb/task_arena.h>

# ifndef CGAL_LINKED_WITH_TBB
#  pragma message(" : Warning: CGAL_LINKED_WITH_TBB not defined: EVERYTHING WILL BE SEQUENTIAL.")
# endif

// # define CGAL_TETRAHEDRAL_REMESHING_USE_SORTED_REFINEMENT_QUEUE
// # define CGAL_TETRAHEDRAL_REMESHING_USE_UNSORTED_REFINEMENT_QUEUE
// # define CGAL_CONCURRENT_TETRAHEDRAL_REMESHING_VERBOSE
// # define CGAL_CONCURRENT_TETRAHEDRAL_REMESHING_VERY_VERBOSE

  // ==========================================================================
  // Concurrency config
  // ==========================================================================

  const char* const CONCURRENT_TETRAHEDRAL_REMESHING_CONFIG_FILENAME = "concurrent_tetrahedral_remeshing_config.cfg";

  // For profiling
# define CGAL_CONCURRENT_TETRAHEDRAL_REMESHING_PROFILING
// # define CGAL_DEBUG_FORCE_SEQUENTIAL_TETRAHEDRAL_REMESHING

#include <thread>

// ==========================================================================
// SEQUENTIAL
// ==========================================================================

#else // !CGAL_CONCURRENT_TETRAHEDRAL_REMESHING

// # define CGAL_TETRAHEDRAL_REMESHING_USE_SORTED_REFINEMENT_QUEUE
// # define CGAL_TETRAHEDRAL_REMESHING_USE_UNSORTED_REFINEMENT_QUEUE

#endif // CGAL_CONCURRENT_TETRAHEDRAL_REMESHING