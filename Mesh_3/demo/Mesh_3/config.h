#ifndef CGAL_DEMO_MESH_3_CONFIG_H
#define CGAL_DEMO_MESH_3_CONFIG_H

#include <CGAL/Mesh_3/Concurrent_mesher_config.h>

// #define CGAL_POLYHEDRON_DEMO_NO_NEF
// #define CGAL_POLYHEDRON_DEMO_NO_SURFACE_MESHER
// #define CGAL_POLYHEDRON_DEMO_NO_PARAMETRIZATION

#define CGAL_MESH_3_VERBOSE
//#define CGAL_MESH_3_VERY_VERBOSE
#define CGAL_MESH_3_IO_VERBOSE

#ifndef CGAL_POLYHEDRON_DEMO_NO_PARAMETRIZATION
#  define CGAL_POLYHEDRON_DEMO_USE_PARAMETRIZATION
#endif

#ifndef CGAL_POLYHEDRON_DEMO_NO_NEF
#  define CGAL_POLYHEDRON_DEMO_USE_NEF
#endif

#ifndef CGAL_POLYHEDRON_DEMO_NO_SURFACE_MESHER
#  define CGAL_POLYHEDRON_DEMO_USE_SURFACE_MESHER
#endif

//#define CGAL_MESH_3_DEMO_ACTIVATE_IMPLICIT_FUNCTIONS
//#define CGAL_MESH_3_DEMO_ACTIVATE_SEGMENTED_IMAGES

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
# define CGAL_CONCURRENT_MESH_3_VERBOSE

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

#   define CGAL_MESH_3_WORKSHARING_USES_TASKS
//#   define CGAL_MESH_3_WORKSHARING_USES_PARALLEL_FOR
//#   define CGAL_MESH_3_WORKSHARING_USES_PARALLEL_DO

    
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

#endif // CONCURRENT_MESH_3
  
#define MESH_3_PROFILING

#endif // CGAL_DEMO_MESH_3_CONFIG_H
