// general includes
#include <iostream>
#include <fstream>
#include <cassert>

// example that uses the Filtered_filtered kernel

#include <CGAL/MP_Float.h>
#include <CGAL/Filtered_kernel.h>

// choose the kernel
#include <CGAL/Simple_cartesian.h>

// inexact kernel
typedef CGAL::Simple_cartesian<double> CK;

// exact kernel
typedef CGAL::Simple_cartesian<CGAL::MP_Float> EK;

typedef CGAL::Filtered_kernel<CK,EK>  Kernel;



// typedefs for the traits and the algorithm

#include <CGAL/Apollonius_graph_2.h>
#include <CGAL/Apollonius_graph_euclidean_traits_2.h>

typedef CGAL::Apollonius_graph_euclidean_traits_2<Kernel> Traits;

// with the second template argument being false, we indicate that
// there is no need to store that hidden weighted points
// one case which is indeed not need is if we do only insertions, like
// in the main program below.
typedef CGAL::Apollonius_graph_2<Traits,false> Apollonius_graph;



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

  return 0;
}


