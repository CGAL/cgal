// standard includes
#include <iostream>
#include <fstream>
#include <cassert>

#include <CGAL/basic.h>

// example that uses the Filtered_exact number type

#include <CGAL/MP_Float.h>
#include <CGAL/Filtered_exact.h>

typedef CGAL::Filtered_exact<double,CGAL::MP_Float> NT;

// choose the kernel
#include <CGAL/Simple_cartesian.h>

typedef CGAL::Simple_cartesian<NT>  Kernel;

// typedefs for the traits and the algorithm

#include <CGAL/Apollonius_graph_2.h>
#include <CGAL/Apollonius_graph_euclidean_traits_2.h>

typedef CGAL::Apollonius_graph_euclidean_traits_2<Kernel> Traits;
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


