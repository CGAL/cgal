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

typedef
CGAL::Segment_Voronoi_diagram_hierarchy_2<Gt>  Segment_Voronoi_diagram;



int main()
{
  std::ifstream ifs("data/sites.cin");
  assert( ifs );

  Segment_Voronoi_diagram          svd;
  Segment_Voronoi_diagram::Site_2  site;

  // read the sites and insert them in the Apollonius graph
  while ( ifs >> site ) {
    if ( site.is_point() ) {
      svd.insert(site.point());
    } else {
      svd.insert(site.source(), site.target());
    }
  }

  // validate the Apollonius graph
  assert( svd.is_valid(true, 1) );
  std::cout << std::endl;

  return 0;
}


