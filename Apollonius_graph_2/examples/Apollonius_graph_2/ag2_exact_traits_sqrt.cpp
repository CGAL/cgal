#include <CGAL/basic.h>

// standard includes
#include <iostream>
#include <fstream>
#include <cassert>


#if defined CGAL_USE_LEDA
#  include <CGAL/leda_real.h>
#elif defined CGAL_USE_CORE
#  include <CGAL/CORE_Expr.h>
#endif


#if defined CGAL_USE_LEDA
// If LEDA is present use leda_real as the exact number type
typedef leda_real NT;

#elif defined CGAL_USE_CORE
// Otherwise if CORE is present use CORE's Expr as the exact number type
typedef CORE::Expr NT;

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
