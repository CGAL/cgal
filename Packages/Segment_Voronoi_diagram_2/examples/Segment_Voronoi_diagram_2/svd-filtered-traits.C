// file: svd-filtered-traits.C
#include <CGAL/basic.h>

// standard includes
#include <iostream>
#include <fstream>
#include <cassert>

// example that uses the filtered traits and
// the segment Voronoi diagram hierarchy

// choose the kernel
#include <CGAL/Simple_cartesian.h>

struct Rep : public CGAL::Simple_cartesian<double> {};

// typedefs for the traits and the algorithm
#include <CGAL/Segment_Voronoi_diagram_hierarchy_2.h>
#include <CGAL/Segment_Voronoi_diagram_filtered_traits_2.h>

struct Gt
  : public CGAL::Segment_Voronoi_diagram_filtered_traits_2<Rep> {};

typedef CGAL::Segment_Voronoi_diagram_hierarchy_2<Gt>  SVD2;


int main()
{
  std::ifstream ifs("data/sites.cin");
  assert( ifs );

  SVD2          svd;
  SVD2::Site_2  site;

  // read the sites and insert them in the segment Voronoi diagram
  while ( ifs >> site ) {
    svd.insert(site);
  }

  // validate the segment Voronoi diagram
  assert( svd.is_valid(true, 1) );

  return 0;
}


