#ifndef CGAL_DEMO_MESH_3_CONFIG_H
#define CGAL_DEMO_MESH_3_CONFIG_H

// #define CGAL_POLYHEDRON_DEMO_NO_NEF
// #define CGAL_POLYHEDRON_DEMO_NO_SURFACE_MESHER
// #define CGAL_POLYHEDRON_DEMO_NO_PARAMETRIZATION

//#define CHECK_AND_DISPLAY_THE_NUMBER_OF_BAD_ELEMENTS_IN_THE_END

//#define CGAL_MESH_3_VERBOSE
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

#define CGAL_MESH_3_DEMO_ACTIVATE_IMPLICIT_FUNCTIONS
//#define CGAL_MESH_3_DEMO_ACTIVATE_SEGMENTED_IMAGES

// ==========================================================================
// MESH_3 GENERAL PARAMETERS
// ==========================================================================

#define CGAL_MESH_3_INITIAL_POINTS_NO_RANDOM_SHOOTING
//#define CGAL_MESHES_DEBUG_REFINEMENT_POINTS

// ==========================================================================
// CONCURRENCY
// ==========================================================================

#ifdef CONCURRENT_MESH_3

# ifndef CGAL_LINKED_WITH_TBB
#   pragma error("CGAL_LINKED_WITH_TBB not defined.")
# endif

# define CGAL_MESH_3_USE_LAZY_SORTED_REFINEMENT_QUEUE
//# define CGAL_MESH_3_USE_LAZY_UNSORTED_REFINEMENT_QUEUE // CJTODO: TEST performance avec et sans
# include <CGAL/Mesh_3/Concurrent_mesher_config.h>

  // ==========================================================================
  // Concurrency activation
  // ==========================================================================

  // In case some code uses CGAL_PROFILE, it needs to be concurrent
# define CGAL_CONCURRENT_PROFILE
# define CGAL_CONCURRENT_MESH_3_VERBOSE
//#define CGAL_CONCURRENT_MESH_3_VERY_VERBOSE

  // ==========================================================================
  // Concurrency config
  // ==========================================================================

  const char * const CONFIG_FILENAME = 
    "D:/INRIA/CGAL/workingcopy/Mesh_3/demo/Mesh_3/concurrent_mesher_config.cfg";

# define CGAL_MESH_3_ADD_OUTSIDE_POINTS_ON_A_FAR_SPHERE
//# define CGAL_MESH_3_ACTIVATE_GRID_INDEX_CACHE_IN_VERTEX // DOES NOT WORK YET
    
  // =================
  // Locking strategy
  // =================
    
# define CGAL_MESH_3_LOCKING_STRATEGY_SIMPLE_GRID_LOCKING

//# define CGAL_MESH_3_CONCURRENT_REFINEMENT_LOCK_ADJ_CELLS // USELESS, FOR TESTS ONLY

  // =====================
  // Worksharing strategy
  // =====================
      
//# define CGAL_MESH_3_WORKSHARING_USES_PARALLEL_FOR
//# define CGAL_MESH_3_WORKSHARING_USES_PARALLEL_DO
# define CGAL_MESH_3_WORKSHARING_USES_TASK_SCHEDULER
# ifdef CGAL_MESH_3_WORKSHARING_USES_TASK_SCHEDULER
//#   define CGAL_MESH_3_TASK_SCHEDULER_WITH_LOCALIZATION_IDS
//#   ifdef CGAL_MESH_3_TASK_SCHEDULER_WITH_LOCALIZATION_IDS
//#     define CGAL_MESH_3_TASK_SCHEDULER_WITH_LOCALIZATION_IDS_SHARED // optional
//#   endif
//#   define CGAL_MESH_3_LOAD_BASED_WORKSHARING // Not recommended
//#   define CGAL_MESH_3_TASK_SCHEDULER_SORTED_BATCHES_WITH_MULTISET
#   define CGAL_MESH_3_TASK_SCHEDULER_SORTED_BATCHES_WITH_SORT
# endif


//#   define CGAL_MESH_3_WORKSHARING_USES_PARALLEL_FOR
//#   define CGAL_MESH_3_WORKSHARING_USES_PARALLEL_DO

  // ==========================================================================
  // Profiling
  // ==========================================================================

  // For abortion profiling, etc.
# define CGAL_CONCURRENT_MESH_3_PROFILING
  
  // ==========================================================================
  // TBB
  // ==========================================================================

  // Use TBB malloc proxy (for all new/delete/malloc/free calls)
  // Highly recommended
# include <tbb/tbbmalloc_proxy.h>



#else // !CONCURRENT_MESH_3

//# define CGAL_MESH_3_USE_LAZY_SORTED_REFINEMENT_QUEUE
# define CGAL_MESH_3_USE_LAZY_UNSORTED_REFINEMENT_QUEUE
# define CGAL_MESH_3_IF_UNSORTED_QUEUE_JUST_SORT_AFTER_SCAN

// For better performance on meshes like fandisk
# define CGAL_MESH_3_ADD_OUTSIDE_POINTS_ON_A_FAR_SPHERE 

#endif // CONCURRENT_MESH_3
  
#define MESH_3_PROFILING

#endif // CGAL_DEMO_MESH_3_CONFIG_H
