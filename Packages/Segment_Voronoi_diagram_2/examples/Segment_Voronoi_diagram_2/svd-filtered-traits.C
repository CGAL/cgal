// file: svd-filtered-traits.C
#include <CGAL/basic.h>

// standard includes
#include <iostream>
#include <fstream>
#include <cassert>

// Workaround for buggy compilers.
#if 0
#ifdef CGAL_CFG_MATCHING_BUG_2
#  define CGAL_IA_CT double
#  define CGAL_IA_PROTECTED true
#  define CGAL_IA_CACHE No_Filter_Cache
#  define CGAL_IA_ET CGAL::MP_Float
#endif
#endif

// example that uses the filtered traits and
// the segment Voronoi diagram hierarchy

// choose the kernel
#include <CGAL/Simple_cartesian.h>

struct Rep : public CGAL::Simple_cartesian<double> {};

// typedefs for the traits and the algorithm

#include <CGAL/Segment_Voronoi_diagram_hierarchy_2.h>
#include <CGAL/Segment_Voronoi_diagram_traits_2.h>

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
    svd.insert(site);
  }

  // validate the Apollonius graph
  assert( svd.is_valid(true, 1) );
  std::cout << std::endl;

  return 0;
}


