#include <CGAL/Random.h>

// example that shows how to add info to input sites and how this is
// propagated using the storage traits with info
//
// the input sites are considered to be colored either red of blue;
// points on the plane belonging to sites of different colors become
// purple.

// standard includes
#include <iostream>
#include <fstream>
#include <cassert>

// an enum representing the color
enum Red_blue {
  RED = 1,
  BLUE = 2,
  PURPLE = 3
};

// output operator for the color
std::ostream&
operator<<(std::ostream& os, const Red_blue& rb)
{
  if ( rb == RED ) { os << "Red"; }
  else if ( rb == BLUE ) { os << "Blue"; }
  else if ( rb == PURPLE ) { os << "Purple"; }
  return os;
}

// functor that defines how to convert color info when:
// 1. constructing the storage site of an endpoint of a segment
// 2. a segment site is split into two sub-segments
struct Red_blue_convert_info
{
  typedef Red_blue      Info;
  typedef const Info&   result_type;

  inline
  const Info& operator()(const Info& info0, bool) const {
    // just return the info of the supporting segment
    return info0;
  }

  inline
  const Info& operator()(const Info& info0, const Info& , bool) const {
    // just return the info of the supporting segment
    return info0;
  }
};


// functor that defines how to merge color info when a site (either
// point or segment) corresponds to point(s) on plane belonging to
// more than one input site
struct Red_blue_merge_info
{
  typedef Red_blue   Info;
  typedef Info       result_type;

  inline
  Info operator()(const Info& info0, const Info& info1) const {
    // if the two sites defining the new site have the same info, keep
    // this common info
    if ( info0 == info1 ) { return info0; }
    // otherwise the new site should be purple
    return PURPLE;
  }
};


// choose the kernel
#include <CGAL/Simple_cartesian.h>

struct Rep : public CGAL::Simple_cartesian<double> {};

// typedefs for the geometric traits, storage traits and the algorithm
#include <CGAL/Segment_Delaunay_graph_2.h>
#include <CGAL/Segment_Delaunay_graph_filtered_traits_2.h>
#include <CGAL/Segment_Delaunay_graph_storage_traits_with_info_2.h>


typedef CGAL::Segment_Delaunay_graph_filtered_traits_2<Rep> Gt;

// define the storage traits with info
typedef
CGAL::Segment_Delaunay_graph_storage_traits_with_info_2<Gt,
							Red_blue,
							Red_blue_convert_info,
							Red_blue_merge_info>
ST;

typedef CGAL::Segment_Delaunay_graph_2<Gt,ST>  SDG2;

typedef SDG2::Finite_vertices_iterator  FVIT;
typedef SDG2::Site_2                    Site_2;


int main()
{

  std::ifstream ifs("data/sitesxx.rb.cin");
  assert( ifs );

  SDG2 sdg;
  Site_2 site;

  // read the sites and their info and insert them in the
  // segment Delaunay graph; print them as you read them
  std::cout << std::endl;
  std::cout << "Input sites:" << std::endl;
  std::cout << "------------" << std::endl;
  while ( ifs >> site ) {
    char c;
    ifs >> c;
    Red_blue info = (c == 'r') ? RED : BLUE;
    std::cout << site << std::flush;
    std::cout << "\r\t\t\t" << info << std::endl;
    sdg.insert(site, info);
  }
  std::cout << std::endl;

  // validate the segment Delaunay graph
  assert( sdg.is_valid(true, 1) );

  // print the sites of the segment Delaunay graph and their info
  std::cout << std::endl << std::endl;
  std::cout << "Output sites:" << std::endl;
  std::cout << "-------------" << std::endl;
  for (FVIT it = sdg.finite_vertices_begin();
       it != sdg.finite_vertices_end(); ++it) {
    std::cout << it->site() << std::flush;
    std::cout << "\r\t\t\t" << it->storage_site().info() << std::endl;
  }
  std::cout << std::endl;

  return 0;
}
