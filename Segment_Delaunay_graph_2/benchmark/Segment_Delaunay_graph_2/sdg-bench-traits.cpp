#include <iostream>
#include <fstream>
#include <cassert>

// define the input kernel
#include <CGAL/Simple_cartesian.h>
typedef CGAL::Simple_cartesian<double>     CK;

// typedefs for the traits and the algorithm
#include <CGAL/Segment_Delaunay_graph_filtered_traits_2.h>
#include <CGAL/Segment_Delaunay_graph_2.h>
#include <CGAL/Segment_Delaunay_graph_hierarchy_2.h>
#include <CGAL/Timer.h>

typedef CGAL::Segment_Delaunay_graph_filtered_traits_2<
          CK,CGAL::Field_with_sqrt_tag>                   Gt_i;
typedef CGAL::Segment_Delaunay_graph_filtered_traits_without_intersections_2<
          CK,CGAL::Field_with_sqrt_tag>                   Gt_wi;

typedef CGAL::Segment_Delaunay_graph_2<Gt_i>              SDG2_i;
typedef CGAL::Segment_Delaunay_graph_2<Gt_wi>             SDG2_wi;
typedef CGAL::Segment_Delaunay_graph_hierarchy_2<Gt_i>    SDG2_Hi;
typedef CGAL::Segment_Delaunay_graph_hierarchy_2<Gt_wi>   SDG2_Hwi;


template <class SDG>
void run_sp(
  const std::vector<CK::Point_2>& points,
  const std::vector< std::pair<std::size_t,std::size_t> >& indices,
  std::string type ){
  CGAL::Timer time;
  time.start();
  SDG sdg;
  sdg.insert_segments( points.begin(),
                       points.end(),
                       indices.begin(),
                       indices.end() );
  time.stop();
  std::cout << type << " with spatial sorting " << time.time() << " "
            << sdg.number_of_input_sites() << std::endl;
}

template <class SDG>
void run_seq(
  const std::vector<CK::Point_2>& points,
  const std::vector< std::pair<std::size_t,std::size_t> >& indices,
  std::string type ){
  CGAL::Timer time;
  time.start();
  SDG sdg;
  std::size_t nbseg = indices.size();
  for (std::size_t i=0; i< nbseg; ++i)
  {
    sdg.insert(
      SDG::Site_2::construct_site_2(
        points[ indices[i].first ],
        points[ indices[i].second ]
      ) );
  }
  time.stop();
  std::cout << type << " sequential " << time.time() << " "
            << sdg.number_of_input_sites() << std::endl;
}
template <class SDG>
void run( const std::vector<CK::Point_2>& points,
          const std::vector< std::pair<std::size_t,std::size_t> >& indices,
          std::string type )
{
  run_seq<SDG>(points, indices, type);
  run_sp<SDG>(points, indices, type);
}

int main(int argc, char** argv) {
  std::ifstream ifs(argv[1]);
  assert( ifs );

  SDG2_wi::Site_2  site;

  std::vector<CK::Point_2> points;
  std::vector<std::pair<std::size_t,std::size_t> > indices;

  //read the number of sites
  std::size_t n;
  ifs >> n;

  //read a close polygon given as its segments
  ifs >> site;
  points.push_back( site.source_of_supporting_site() );

  std::size_t k=0;
  while ( --n != 0) {
    ifs >> site;
    points.push_back( site.source_of_supporting_site() );
    indices.push_back( std::make_pair(k, k+1) );
    ++k;
  }
  indices.push_back( std::make_pair(k, 0) );
  ifs.close();

  std::cout << "read a " << points.size() << " vertices polygon\n";

  run<SDG2_i>(points, indices, "SDG2_i");
  run<SDG2_wi>(points, indices, "SDG2_wi");
  run<SDG2_Hi>(points, indices, "SDG2_Hi");
  run<SDG2_Hwi>(points, indices, "SDG2_Hwi");

  return 0;
}
