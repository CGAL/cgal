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

// SAME COMMENT HERE; I WILL RESTORE IT WHEN FILTERED KERNEL WORKS
//typedef CGAL::Filtered_kernel<CK,EK>  Kernel;
typedef CGAL::Filtered_kernel<CK>  Kernel;



// typedefs for the traits and the algorithm

#include <CGAL/Apollonius_graph_hierarchy_2.h>
#include <CGAL/Apollonius_graph_euclidean_traits_2.h>

typedef CGAL::Apollonius_graph_euclidean_traits_2<Kernel> Traits;

// now we use the Apollonius graph hierarchy.
// the hierarchy is faster for inputs consisting of about more than
// 100,000 weighted points
typedef CGAL::Apollonius_graph_hierarchy_2<Traits> Apollonius_graph;



int main()
{
  std::ifstream ifs("data/hierarchy.cin");
  assert( ifs );

  Apollonius_graph ag;
  Apollonius_graph::Site_2 site;

  // read the weighted points and insert them in the Apollonius graph
  while ( ifs >> site ) {
    ag.insert(site);
  }

  // validate the Apollonius graph
  assert( ag.is_valid(true, 1) );

  return 0;
}


