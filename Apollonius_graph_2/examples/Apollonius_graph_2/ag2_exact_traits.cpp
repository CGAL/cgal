// standard includes
#include <iostream>
#include <fstream>
#include <cassert>

// the number type
#include <CGAL/MP_Float.h>


// example that uses an exact number type

typedef CGAL::MP_Float NT;

// *** WARNING ***
// The use of a kernel based on an exact number type is highly inefficient.
// It is used in this example primarily for illustration purposes.
// In an efficiency critical context, and/or for the purposes of
// benchmarking the Apollonius_graph_filtered_traits_2<> class should
// be used.

// choose the kernel
#include <CGAL/Simple_cartesian.h>

typedef CGAL::Simple_cartesian<NT>  Kernel;

// typedefs for the traits and the algorithm

#include <CGAL/Apollonius_graph_2.h>
#include <CGAL/Apollonius_graph_traits_2.h>

typedef CGAL::Apollonius_graph_traits_2<Kernel>   Traits;
typedef CGAL::Apollonius_graph_2<Traits>          Apollonius_graph;



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
