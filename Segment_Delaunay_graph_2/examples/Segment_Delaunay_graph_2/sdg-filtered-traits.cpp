// standard includes
#include <iostream>
#include <fstream>
#include <cassert>

// example that uses the filtered traits and
// the segment Delaunay graph hierarchy

// choose the kernel
#include <CGAL/Simple_cartesian.h>

typedef CGAL::Simple_cartesian<double> K;

// typedefs for the traits and the algorithm
#include <CGAL/Segment_Delaunay_graph_hierarchy_2.h>
#include <CGAL/Segment_Delaunay_graph_filtered_traits_2.h>

typedef CGAL::Segment_Delaunay_graph_filtered_traits_2<K> Gt;

typedef CGAL::Segment_Delaunay_graph_hierarchy_2<Gt>  SDG2;


int main()
{
  std::ifstream ifs("data/sites.cin");
  assert( ifs );

  SDG2          sdg;
  SDG2::Site_2  site =
    SDG2::Site_2::construct_site_2(K::Point_2(CGAL::ORIGIN));

  // read the sites and insert them in the segment Delaunay graph
  while ( ifs >> site ) {
    sdg.insert(site);
  }

  // validate the segment Delaunay graph
  assert( sdg.is_valid(true, 1) );

  return 0;
}
