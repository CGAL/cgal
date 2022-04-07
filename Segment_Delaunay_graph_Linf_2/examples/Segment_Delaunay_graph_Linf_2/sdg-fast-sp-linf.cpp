// standard includes
#include <iostream>
#include <fstream>
#include <cassert>

// example that uses the filtered traits,
// the segment Delaunay graph and the spatial sorting

// choose the kernel
#include <CGAL/Simple_cartesian.h>

typedef CGAL::Simple_cartesian<double> K;

// typedefs for the traits and the algorithm
#include <CGAL/Segment_Delaunay_graph_Linf_2.h>
#include <CGAL/Segment_Delaunay_graph_Linf_filtered_traits_2.h>

typedef CGAL::Segment_Delaunay_graph_Linf_filtered_traits_without_intersections_2<K> Gt;

typedef CGAL::Segment_Delaunay_graph_Linf_2<Gt>  SDG2;


int main()
{
  std::ifstream ifs("data/sites.cin");
  assert( ifs );

  SDG2          sdg;
  SDG2::Site_2  site;

  std::vector<SDG2::Site_2> sites;
  // read the sites
  while ( ifs >> site ) {
    sites.push_back(site);
  }

  //insert the sites all at once using spatial sorting to speed the insertion
  sdg.insert( sites.begin(), sites.end(),CGAL::Tag_true() );

  // validate the segment Delaunay graph
  assert( sdg.is_valid(true, 1) );

  return 0;
}
