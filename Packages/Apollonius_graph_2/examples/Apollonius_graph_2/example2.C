// standard includes
#include <iostream>
#include <fstream>
#include <cassert>


#include <CGAL/basic.h>


#if defined CGAL_USE_LEDA
#  include <CGAL/leda_real.h>
#elif defined CGAL_USE_CORE
#  include <CGAL_Expr.h>
#endif

#if defined CGAL_USE_LEDA || defined CGAL_USE_CORE

// Workaround for buggy compilers.
#ifdef CGAL_CFG_MATCHING_BUG_2
#  define CGAL_IA_CT double
#  define CGAL_IA_PROTECTED true
#  define CGAL_IA_CACHE No_Filter_Cache
#  ifdef CGAL_USE_LEDA
#    define CGAL_IA_ET leda_real
#  elif defined CGAL_USE_CORE
#    define CGAL_IA_ET CORE::Expr
#  endif
#endif

#include <CGAL/Filtered_exact.h>


#endif


#if defined CGAL_USE_LEDA
// If LEDA is present use leda_real as the exact number type for
// Filtered_exact
typedef CGAL::Filtered_exact<double,leda_real> NT;

#elif defined CGAL_USE_CORE
// Othwrwise if CORE is present use CORE's Expr as the exact number
// type for Filtered_exact

typedef CGAL::Filtered_exact<double,CORE::Expr> NT;

#else

// Otherwise just use double. This may cause numerical errors but it
// is still worth doing it to show how to define correctly the traits
// class
typedef double NT;

#endif


#include <CGAL/Simple_cartesian.h>

typedef CGAL::Simple_cartesian<NT>  Kernel;


// typedefs for the traits and the algorithm

#include <CGAL/Apollonius_graph_2.h>
#include <CGAL/Apollonius_graph_euclidean_traits_2.h>

// the traits class is now going to assume that the operations
// +,-,*,/ and sqrt are supported exactly
typedef
CGAL::Apollonius_graph_euclidean_traits_2<Kernel,CGAL::Sqrt_field_tag>
Traits;

typedef CGAL::Apollonius_graph_2<Traits> Apollonius_graph;



int main()
{
  std::ifstream ifs("data/sites.cin");
  assert( ifs );

  Apollonius_graph ag;
  Apollonius_graph::Apollonius_site_2 site;

  // read the weighted points and insert them in the Apollonius graph
  while ( ifs >> site ) {
    ag.insert(site);
  }

  // validate the Apollonius graph
  assert( ag.is_valid(true, 1) );
  std::cout << std::endl;

  return 0;
}
