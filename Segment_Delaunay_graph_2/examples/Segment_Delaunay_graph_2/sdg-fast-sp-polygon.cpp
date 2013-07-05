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
#include <CGAL/Segment_Delaunay_graph_2.h>
#include <CGAL/Segment_Delaunay_graph_filtered_traits_2.h>

typedef CGAL::Segment_Delaunay_graph_filtered_traits_without_intersections_2<K> Gt;
typedef CGAL::Segment_Delaunay_graph_2<Gt>  SDG2;


int main()
{
  std::ifstream ifs("data/sites.cin");
  assert( ifs );

  //polygon points
  std::vector<Gt::Point_2> points;
  //segments of the polygon as a pair of point indices
  std::vector<std::pair<std::size_t,std::size_t> > indices;

  SDG2::Site_2 site;
  //read a close polygon given by its segments
  // s x0 y0 x1 y1
  // s x1 y1 x2 y2
  // ...
  // s xn yn x0 y0
  ifs >> site;
  assert( site.is_segment() );
  points.push_back( site.source_of_supporting_site() );

  std::size_t k=0;
  while ( ifs >> site ) {
    assert( site.is_segment() );
    points.push_back( site.source_of_supporting_site() );
    indices.push_back( std::make_pair(k, k+1) );
    ++k;
  }
  indices.push_back( std::make_pair(k, 0) );
  ifs.close();

  SDG2          sdg;

  //insert the polygon segments all at once using spatial sorting to speed the insertion
  sdg.insert_segments( points.begin(), points.end(), indices.begin(), indices.end() );

  // validate the segment Delaunay graph
  assert( sdg.is_valid(true, 1) );

  return 0;
}
