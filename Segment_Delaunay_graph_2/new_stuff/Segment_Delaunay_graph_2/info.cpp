// file: info.cpp
#define WITH_INFO

#include <CGAL/basic.h>

// standard includes
#include <iostream>
#include <fstream>
#include <cassert>

// example that uses the filtered traits and
// the segment Delaunay graph hierarchy

// choose the kernel
#include <CGAL/Simple_cartesian.h>

struct Rep : public CGAL::Simple_cartesian<double> {};

// typedefs for the traits and the algorithm
#include <CGAL/Segment_Delaunay_graph_2.h>
#include <CGAL/Segment_Delaunay_graph_hierarchy_2.h>

#include <CGAL/Segment_Delaunay_graph_filtered_traits_2.h>

typedef
CGAL::Segment_Delaunay_graph_filtered_traits_2<Rep>
Traits_x;

typedef
CGAL::Segment_Delaunay_graph_filtered_traits_without_intersections_2<Rep>
Traits_no_x;

struct Gt : public Traits_x {};
//struct Gt : public Traits_no_x {};

#ifdef WITH_INFO
#include "Red_blue_info.h"
#include <CGAL/Storage_traits_with_info_2.h>

typedef CGAL::Storage_traits_2<Gt>   ST_base;
typedef CGAL::Storage_traits_with_info_2<ST_base,
					 ::Red_blue,
					 ::Red_blue_convert_info,
					 ::Red_blue_merge_info>
ST;
#else
typedef CGAL::Storage_traits_2<Gt> ST;
#endif

typedef CGAL::Tag_true  STag;
//typedef CGAL::Tag_false STag;

typedef CGAL::Segment_Delaunay_graph_2<Gt,ST>  SDG2;
//typedef CGAL::Segment_Delaunay_graph_hierarchy_2<Gt,STag,ST>  SDG2;

template<class SDG>
bool test(SDG& sdg, char* fname, bool read_info = false)
{
  std::ifstream ifs(fname);
  assert( ifs );

  sdg.clear();
  typename SDG::Site_2  site;

  // read the sites and insert them in the segment Delaunay graph
  std::cout << "Input:" << std::endl;
#ifdef WITH_INFO
  Random_red_blue random_red_blue;
#endif
  while ( ifs >> site ) {
#ifdef WITH_INFO
    Red_blue info;
    if ( read_info ) {
      char c;
      ifs >> c;
      info = (c == 'r') ? RED : BLUE;
    } else {
      info = *random_red_blue;
    }
    std::cout << "SITE TO BE INSERTED: "
	      << site << " " << info << std::endl;
    sdg.insert(site, info);
#else
    std::cout << "SITE TO BE INSERTED: "
	      << site << std::endl;
    sdg.insert(site);
#endif
  }
  std::cout << std::endl;

  typedef typename SDG::Finite_vertices_iterator FVIT;
  for (FVIT it = sdg.finite_vertices_begin();
       it != sdg.finite_vertices_end(); ++it) {
#ifdef WITH_INFO
    std::cout << it->site() << " "
	      << it->storage_site().info() << std::endl;
#else
    std::cout << it->site() << std::endl;
#endif
  }

  // validate the segment Delaunay graph
  return sdg.is_valid(true, 1);
}

void print_separator()
{
  std::cout << std::endl << std::endl;
  std::cout << "=============================";
  std::cout << "=============================";
  std::cout << std::endl;
  std::cout << "=============================";
  std::cout << "=============================";
  std::cout << std::endl << std::endl;
}

int main()
{
  bool b;
  SDG2 sdg;

  print_separator();
  b = test(sdg, "data/sites.cin");
  assert( b );

  print_separator();
  b = test(sdg, "data/sitesx2.cin");
  assert( b );

  print_separator();
  b = test(sdg, "data/bizarre.cin");
  assert( b );  

  print_separator();
  b = test(sdg, "data/sitesx.rb.cin", true);
  assert( b );

  print_separator();
  b = test(sdg, "data/sitesxx.rb.cin", true);
  assert( b );

  print_separator();

  return 0;
}


