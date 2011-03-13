// standard includes
#include <iostream>
#include <fstream>
#include <cassert>

// example that uses the filtered traits

#include <CGAL/MP_Float.h>
#include <CGAL/Simple_cartesian.h>

// constructions kernel (inexact)
typedef CGAL::Simple_cartesian<double> CK;

// exact kernel
typedef CGAL::Simple_cartesian<CGAL::MP_Float> EK;


// typedefs for the traits and the algorithm

#include <CGAL/Apollonius_graph_hierarchy_2.h>
#include <CGAL/Apollonius_graph_filtered_traits_2.h>


// Type definition for the traits class.
// In this example we explicitly define the exact kernel. We also
// explicitly define what operations to use for the evaluation of the
// predicates and constructions, when the filtering and the exact
// kernels are used respectively.
// Note that the operations allowed for the filtering and the
// constructions (field operations plus square roots) are different
// from the operations allowed when the exact kernel is used (ring
// operations).
typedef CGAL::Field_with_sqrt_tag  CM;
typedef CGAL::Integral_domain_without_division_tag        EM;
typedef CGAL::Apollonius_graph_filtered_traits_2<CK,CM,EK,EM> Traits;

// Now we use the Apollonius graph hierarchy.
// The hierarchy is faster for inputs consisting of about more than
// 1,000 sites
typedef CGAL::Apollonius_graph_hierarchy_2<Traits> Apollonius_graph;



int main()
{
  std::ifstream ifs("data/hierarchy.cin");
  assert( ifs );

  Apollonius_graph ag;
  Apollonius_graph::Site_2 site;

  // read the sites and insert them in the Apollonius graph
  while ( ifs >> site ) {
    ag.insert(site);
  }

  // validate the Apollonius graph
  assert( ag.is_valid(true, 1) );

  return 0;
}
