// general includes
#include <iostream>
#include <fstream>
#include <cassert>

#include <CGAL/Filtered_exact.h>


#ifdef CGAL_USE_LEDA

// If LEDA is present use leda_real as the exact number type for
// Filtered_exact
#include <CGAL/leda_real.h>
typedef CGAL::Filtered_exact<double,leda_real> NT;

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



int main(int argc, char* argv[])
{
  assert( argc >= 2 );

  std::ifstream ifs(argv[1]);
  assert( ifs );

  Apollonius_graph ag;
  Apollonius_graph::Weighted_point wp;

  // read the weighted points and insert them in the Apollonius graph
  while ( ifs >> wp ) {
    ag.insert(wp);
  }

  // validate the Apollonius graph
  assert( ag.is_valid(true, 1) );
  std::cout << std::endl;

  return 0;
}
