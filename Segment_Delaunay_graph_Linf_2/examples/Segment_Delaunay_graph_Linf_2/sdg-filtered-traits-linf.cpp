// standard includes
#include <iostream>
#include <fstream>
#include <cassert>

// example that uses the filtered traits and
// the segment Delaunay graph hierarchy

// choose the kernel
#include <CGAL/Simple_cartesian.h>

typedef CGAL::Simple_cartesian<double> Rep;

// typedefs for the traits and the algorithm
#include <CGAL/Segment_Delaunay_graph_Linf_hierarchy_2.h>
#include <CGAL/Segment_Delaunay_graph_Linf_filtered_traits_2.h>

typedef CGAL::Segment_Delaunay_graph_Linf_filtered_traits_2<Rep> Gt;

typedef CGAL::Segment_Delaunay_graph_Linf_hierarchy_2<Gt>  SDG2;


int main( int argc, char *argv[] )
{
  if ( ! (( argc == 1 ) || (argc == 2)) ) {
    std::cout <<"usage: "<< argv[0] <<" [filename]\n";
  }

  std::ifstream ifs( (argc == 1) ? "data/sites.cin" : argv[1] );
  assert( ifs );

  SDG2          sdg;
  SDG2::Site_2  site =
    SDG2::Site_2::construct_site_2(Rep::Point_2(CGAL::ORIGIN));

  // read the sites and insert them in the segment Delaunay graph
  while ( ifs >> site ) {
    sdg.insert(site);
  }

  // validate the segment Delaunay graph
  assert( sdg.is_valid(true, 1) );

  return 0;
}
