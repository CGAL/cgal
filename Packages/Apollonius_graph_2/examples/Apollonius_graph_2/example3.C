// standard includes
#include <iostream>
#include <fstream>
#include <cassert>

#include <CGAL/basic.h>

// example that uses the Filtered_kernel

#include <CGAL/MP_Float.h>
#include <CGAL/Filtered_kernel.h>

// choose the kernel
#include <CGAL/Simple_cartesian.h>

// inexact kernel
typedef CGAL::Simple_cartesian<double> CK;

// exact kernel
typedef CGAL::Simple_cartesian<CGAL::MP_Float> EK;


typedef CGAL::Filtered_kernel<CK>  Kernel;


// typedefs for the traits and the algorithm

#include <CGAL/Apollonius_graph_2.h>
#include <CGAL/Apollonius_graph_euclidean_traits_2.h>

typedef CGAL::Apollonius_graph_euclidean_traits_2<Kernel> Traits;

// with the second template argument being false, we indicate that
// there is no need to store that hidden weighted points;
// one case where this is indeed not needed is when we only do
// insertions, like in the main program below.
typedef CGAL::Apollonius_graph_2<Traits,false> Apollonius_graph;



int main(int argc, char* argv[])
{
  assert( argc >= 2 );

  std::ifstream ifs(argv[1]);
  assert( ifs );

  Apollonius_graph ag;
  Apollonius_graph::Weighted_point_2 wp;

  // read the weighted points and insert them in the Apollonius graph
  while ( ifs >> wp ) {
    ag.insert(wp);
  }

  // validate the Apollonius graph
  assert( ag.is_valid(true, 1) );
  std::cout << std::endl;

  return 0;
}


