#include <CGAL/config.h>

// standard includes
#include <iostream>
#include <fstream>
#include <cassert>

#if defined(CGAL_USE_CORE) || defined(CGAL_USE_LEDA)
#  include <CGAL/Exact_algebraic.h>
#endif

// *** WARNING ***
// The use of a kernel based on an exact number type is highly inefficient.
// It is used in this example primarily for illustration purposes.
// In an efficiency critical context, and/or for the purposes of
// benchmarking the Apollonius_graph_filtered_traits_2<> class should
// be used.

#if defined(CGAL_USE_CORE) || defined(CGAL_USE_LEDA)
// If LEDA is present use leda_real as the exact number type
// Otherwise if CORE is present use CORE's Expr as the exact number type
typedef CGAL::Exact_algebraic NT;

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
#include <CGAL/Apollonius_graph_traits_2.h>

// the traits class is now going to assume that the operations
// +,-,*,/ and sqrt are supported exactly
typedef
CGAL::Apollonius_graph_traits_2<Kernel,CGAL::Field_with_sqrt_tag> Traits;

typedef CGAL::Apollonius_graph_2<Traits> Apollonius_graph;



int main()
{
  std::ifstream ifs("data/sites.cin");
  assert( ifs );

  Apollonius_graph ag;
  Apollonius_graph::Site_2 site;

  // read the sites and insert them in the Apollonius graph
  while ( ifs >> site ) {
    ag.insert(site);
  }

  // validate the Apollonius graph
  assert( ag.is_valid(true, 1) );
  std::cout << std::endl;

  return 0;
}
