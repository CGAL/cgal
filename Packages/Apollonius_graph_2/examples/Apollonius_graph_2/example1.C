// general includes
#include <iostream>
#include <fstream>
#include <cassert>

// example that uses the Filtered_exact number type

#include <CGAL/MP_Float.h>
#include <CGAL/Filtered_exact.h>

// I SHOULD BE ABLE TO DO THE FOLLOWING INSTEAD OF JUST PUTTING ONLY
//ONE ARGUMENT
typedef CGAL::Filtered_exact<double,CGAL::MP_Float> NT;

// choose the kernel
#include <CGAL/Simple_cartesian.h>

typedef CGAL::Simple_cartesian<NT>  Kernel;

// typedefs for the traits and the algorithm

#include <CGAL/Apollonius_graph_2.h>
#include <CGAL/Apollonius_graph_euclidean_traits_2.h>

typedef CGAL::Apollonius_graph_euclidean_traits_2<Kernel> Traits;
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


