#ifndef CGAL_DEMO_MESH_3_CONFIG_H
#define CGAL_DEMO_MESH_3_CONFIG_H

#define  BOOST_PARAMETER_MAX_ARITY 12

// CGAL_MESH_3_PROTECTION_DEBUG:
// -  1 : display debug messages
// -  2 : dump file `polylines_graph.polylines.txt` and
//        `edges-graph.polylines.txt` from polyhedral and image domains
// -  4 : dump c3t3 in case of a bug
// -  8 : dump c3t3 at various stages of the protection
// - 16 : more precise debug messages
//#define CGAL_MESH_3_PROTECTION_DEBUG 255


//#define CGAL_PROFILE

// #define CGAL_POLYHEDRON_DEMO_NO_NEF
// #define CGAL_POLYHEDRON_DEMO_NO_SURFACE_MESHER
// #define CGAL_POLYHEDRON_DEMO_NO_PARAMETRIZATION

//#define CGAL_MESH_3_VERBOSE
//#define CGAL_MESH_3_PERTURBER_HIGH_VERBOSITY
//#define CGAL_MESH_3_EXUDER_VERBOSE
//#define CGAL_MESH_3_EXUDER_HIGH_VERBOSITY
//#define CGAL_MESH_3_VERY_VERBOSE
//#define CGAL_MESH_3_IO_VERBOSE

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

//#define CGAL_MESH_3_DEMO_BIGGER_HISTOGRAM_WITH_WHITE_BACKGROUNG

// If you define this, implicit function and segmented images won't be available
//#define CGAL_MESH_3_DEMO_ACTIVATE_SHARP_FEATURES_IN_POLYHEDRAL_DOMAIN
#ifndef CGAL_MESH_3_DEMO_ACTIVATE_SHARP_FEATURES_IN_POLYHEDRAL_DOMAIN
# define CGAL_MESH_3_DEMO_ACTIVATE_IMPLICIT_FUNCTIONS
# define CGAL_MESH_3_DEMO_ACTIVATE_SEGMENTED_IMAGES
#endif

//#define CGAL_MESH_3_DEMO_DONT_COUNT_TETS_ADJACENT_TO_SHARP_FEATURES_FOR_HISTOGRAM

// Optimizers
//#define CGAL_MESH_3_DEMO_DISABLE_ODT
//#define CGAL_MESH_3_DEMO_DISABLE_LLOYD
//#define CGAL_MESH_3_DEMO_DISABLE_PERTURBER
//#define CGAL_MESH_3_DEMO_DISABLE_EXUDER

// ==========================================================================
// MESH_3 GENERAL PARAMETERS
// ==========================================================================

//#define CGAL_MESH_3_USE_OLD_SURFACE_RESTRICTED_DELAUNAY_UPDATE // WARNING: VERY SLOW
//#define CGAL_MESH_3_INITIAL_POINTS_NO_RANDOM_SHOOTING
//#define CGAL_MESHES_DEBUG_REFINEMENT_POINTS
//#define CHECK_AND_DISPLAY_THE_NUMBER_OF_BAD_ELEMENTS_IN_THE_END

// ==========================================================================
// ==========================================================================
// CONCURRENT MESH_3?
// ==========================================================================
// ==========================================================================

#ifdef CGAL_CONCURRENT_MESH_3

# ifndef CGAL_LINKED_WITH_TBB
#   pragma message(" : Warning: CGAL_LINKED_WITH_TBB not defined: EVERYTHING WILL BE SEQUENTIAL.")
# endif

# define CGAL_MESH_3_USE_LAZY_SORTED_REFINEMENT_QUEUE // default behavior
//# define CGAL_MESH_3_USE_LAZY_UNSORTED_REFINEMENT_QUEUE
# define CGAL_MESH_3_IF_UNSORTED_QUEUE_JUST_SORT_AFTER_SCAN
# include <CGAL/Mesh_3/Concurrent_mesher_config.h>

  // ==========================================================================
  // Verbose
  // ==========================================================================

//# define CGAL_CONCURRENT_MESH_3_VERBOSE
//#define CGAL_CONCURRENT_MESH_3_VERY_VERBOSE

  // ==========================================================================
  // Concurrency config
  // ==========================================================================

  const char * const CONFIG_FILENAME = "concurrent_mesher_config.cfg";

  // =====================
  // Worksharing strategy
  // =====================

//# define CGAL_MESH_3_LOAD_BASED_WORKSHARING // Not recommended

  // ==========================================================================
  // Profiling
  // ==========================================================================

  // For abortion profiling, etc.
//# define CGAL_CONCURRENT_MESH_3_PROFILING
  // Debugging
//# define CGAL_DEBUG_FORCE_SEQUENTIAL_MESH_REFINEMENT


// ==========================================================================
// ==========================================================================
// SEQUENTIAL MESH_3?
// ==========================================================================
// ==========================================================================

#else // !CGAL_CONCURRENT_MESH_3

//# define CGAL_MESH_3_USE_LAZY_SORTED_REFINEMENT_QUEUE
//# define CGAL_MESH_3_USE_LAZY_UNSORTED_REFINEMENT_QUEUE
# define CGAL_MESH_3_IF_UNSORTED_QUEUE_JUST_SORT_AFTER_SCAN

#endif // CGAL_CONCURRENT_MESH_3

//#define CGAL_MESH_3_PROFILING

#endif // CGAL_DEMO_MESH_3_CONFIG_H
