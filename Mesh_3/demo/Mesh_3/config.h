#ifndef CGAL_DEMO_MESH_3_CONFIG_H
#define CGAL_DEMO_MESH_3_CONFIG_H

// #define CGAL_POLYHEDRON_DEMO_NO_NEF
// #define CGAL_POLYHEDRON_DEMO_NO_SURFACE_MESHER
// #define CGAL_POLYHEDRON_DEMO_NO_PARAMETRIZATION

//#define CGAL_MESH_3_VERBOSE
//#define CGAL_MESH_3_PERTURBER_HIGH_VERBOSITY
//#define CGAL_MESH_3_EXUDER_VERBOSE
//#define CGAL_MESH_3_EXUDER_HIGH_VERBOSITY
//#define CGAL_MESH_3_VERY_VERBOSE
#define CGAL_MESH_3_IO_VERBOSE
#define CGAL_DEBUG_GLOBAL_LOCK_DS // CJTODO TEMP

//#define SHOW_REMAINING_BAD_ELEMENT_IN_RED

#ifndef CGAL_POLYHEDRON_DEMO_NO_PARAMETRIZATION
#  define CGAL_POLYHEDRON_DEMO_USE_PARAMETRIZATION
#endif

#ifndef CGAL_POLYHEDRON_DEMO_NO_NEF
#  define CGAL_POLYHEDRON_DEMO_USE_NEF
#endif

#ifndef CGAL_POLYHEDRON_DEMO_NO_SURFACE_MESHER
#  define CGAL_POLYHEDRON_DEMO_USE_SURFACE_MESHER
#endif

#define CGAL_MESH_3_DEMO_BIGGER_HISTOGRAM_WITH_WHITE_BACKGROUNG

// If you define this, implicit function and segmented images won't be available
//#define CGAL_MESH_3_DEMO_ACTIVATE_SHARP_FEATURES_IN_POLYHEDRAL_DOMAIN
#ifndef CGAL_MESH_3_DEMO_ACTIVATE_SHARP_FEATURES_IN_POLYHEDRAL_DOMAIN
# define CGAL_MESH_3_DEMO_ACTIVATE_IMPLICIT_FUNCTIONS
# define CGAL_MESH_3_DEMO_ACTIVATE_SEGMENTED_IMAGES
#endif

//#define CGAL_MESH_3_DEMO_DONT_COUNT_TETS_ADJACENT_TO_SHARP_FEATURES_FOR_HISTOGRAM

// Debugging
//#define CGAL_DEBUG_FORCE_SEQUENTIAL_MESH_REFINEMENT

// Optimizers
//#define CGAL_MESH_3_DEMO_DISABLE_ODT
//#define CGAL_MESH_3_DEMO_DISABLE_LLOYD
//#define CGAL_MESH_3_DEMO_DISABLE_PERTURBER
//#define CGAL_MESH_3_DEMO_DISABLE_EXUDER

// ==========================================================================
// MESH_3 GENERAL PARAMETERS
// ==========================================================================

//#define CGAL_MESH_3_USE_OLD_SURFACE_RESTRICTED_DELAUNAY_UPDATE // WARNING: VERY SLOW
#define CGAL_MESH_3_INITIAL_POINTS_NO_RANDOM_SHOOTING
#define CGAL_MESH_3_ADD_OUTSIDE_POINTS_ON_A_FAR_SPHERE
//#define CGAL_MESHES_DEBUG_REFINEMENT_POINTS
//#define CHECK_AND_DISPLAY_THE_NUMBER_OF_BAD_ELEMENTS_IN_THE_END

// ==========================================================================
// ==========================================================================
// CONCURRENT MESH_3?
// ==========================================================================
// ==========================================================================

#ifdef CONCURRENT_MESH_3

# ifndef CGAL_LINKED_WITH_TBB
#   pragma message(" : Warning: CGAL_LINKED_WITH_TBB not defined: EVERYTHING WILL BE SEQUENTIAL.")
# endif

# define CGAL_MESH_3_USE_LAZY_SORTED_REFINEMENT_QUEUE // default behavior
//# define CGAL_MESH_3_USE_LAZY_UNSORTED_REFINEMENT_QUEUE
# define CGAL_MESH_3_IF_UNSORTED_QUEUE_JUST_SORT_AFTER_SCAN
# include <CGAL/Mesh_3/Concurrent_mesher_config.h>

  // ==========================================================================
  // Concurrency activation
  // ==========================================================================

# define CGAL_CONCURRENT_MESH_3_VERBOSE
//#define CGAL_CONCURRENT_MESH_3_VERY_VERBOSE

  // ==========================================================================
  // Concurrency config
  // ==========================================================================

  const char * const CONFIG_FILENAME = 
    "D:/INRIA/CGAL/workingcopy/Mesh_3/demo/Mesh_3/concurrent_mesher_config.cfg";
  
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
//#   define CGAL_MESH_3_LOAD_BASED_WORKSHARING // Not recommended
#   define CGAL_MESH_3_TASK_SCHEDULER_SORTED_BATCHES_WITH_MULTISET
//#   define CGAL_MESH_3_TASK_SCHEDULER_SORTED_BATCHES_WITH_SORT
# endif


//#   define CGAL_MESH_3_WORKSHARING_USES_PARALLEL_FOR
//#   define CGAL_MESH_3_WORKSHARING_USES_PARALLEL_DO

  // ==========================================================================
  // Profiling
  // ==========================================================================

  // For abortion profiling, etc.
# define CGAL_CONCURRENT_MESH_3_PROFILING


// ==========================================================================
// ==========================================================================
// SEQUENTIAL MESH_3?
// ==========================================================================
// ==========================================================================

#else // !CONCURRENT_MESH_3

//# define CGAL_MESH_3_USE_LAZY_SORTED_REFINEMENT_QUEUE
//# define CGAL_MESH_3_USE_LAZY_UNSORTED_REFINEMENT_QUEUE
# define CGAL_MESH_3_IF_UNSORTED_QUEUE_JUST_SORT_AFTER_SCAN

#endif // CONCURRENT_MESH_3
  
#define MESH_3_PROFILING

#endif // CGAL_DEMO_MESH_3_CONFIG_H
