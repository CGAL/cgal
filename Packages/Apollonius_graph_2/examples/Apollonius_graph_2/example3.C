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
#include <CGAL/Apollonius_graph_data_structure_2.h>
#include <CGAL/Apollonius_graph_vertex_base_2.h>
#include <CGAL/Apollonius_graph_face_base_2.h>
#include <CGAL/Apollonius_graph_euclidean_traits_2.h>

typedef CGAL::Apollonius_graph_euclidean_traits_2<Kernel> Traits;

// with the second template argument in the vertex base class being
// false, we indicate that there is no need to store that hidden
// weighted points;
// one case where this is indeed not needed is when we only do
// insertions, like in the main program below.
typedef CGAL::Apollonius_graph_vertex_base_2<Traits,false>   Vb;
typedef CGAL::Apollonius_graph_face_base_2<Traits>           Fb;
typedef CGAL::Apollonius_graph_data_structure_2<Vb,Fb>       Agds;
typedef CGAL::Apollonius_graph_2<Traits,Agds>         Apollonius_graph;



int main()
{
  std::ifstream ifs("data/sites.cin");
  assert( ifs );

  Apollonius_graph ag;
  Apollonius_graph::Site_2 site;

  // read the weighted points and insert them in the Apollonius graph
  while ( ifs >> site ) {
    ag.insert(site);
  }

  // validate the Apollonius graph
  assert( ag.is_valid(true, 1) );
  std::cout << std::endl;

  return 0;
}


