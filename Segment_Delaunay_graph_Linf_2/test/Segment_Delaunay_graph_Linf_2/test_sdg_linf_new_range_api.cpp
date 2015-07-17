#ifndef CGAL_SDG_VERBOSE
#define CGAL_SDG_DEBUG(a)
#else
#define CGAL_SDG_DEBUG(a) { a }
#endif

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Segment_Delaunay_graph_Linf_filtered_traits_2.h>
#include <CGAL/Segment_Delaunay_graph_Linf_2.h>
#include <CGAL/Gmpq.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Field_with_sqrt_tag MTag;
typedef CGAL::Integral_domain_without_division_tag EMTag;
typedef CGAL::Simple_cartesian<CGAL::Gmpq> EK;

typedef CGAL::Segment_Delaunay_graph_Linf_filtered_traits_without_intersections_2<K, MTag, EK, EMTag>  Gt;
typedef CGAL::Segment_Delaunay_graph_Linf_2<Gt>  SDG;
typedef SDG::Point_2 Point;
typedef SDG::Site_2 Site_2;


int main()
{
  std::vector<Point> points;
  points.push_back(Point(0,0));
  points.push_back(Point(0,1));
  points.push_back(Point(1,1));
  points.push_back(Point(2,3));
  points.push_back(Point(0,8));

  {
  std::cout <<"Testing insertion of a range of points\n";
  SDG sdg;
  sdg.insert_points(points.begin(), points.end());
  }

  {
  std::cout <<"Testing insertion of a range of segments from points and indices\n";
  SDG sdg;
  std::vector< std::pair<int,int> > indices;
  indices.push_back( std::make_pair(0,1) );
  indices.push_back( std::make_pair(1,2) );
  indices.push_back( std::make_pair(2,3) );
  indices.push_back( std::make_pair(3,4) );
  indices.push_back( std::make_pair(4,0) );
  sdg.insert_segments( points.begin(), points.end(),
                       indices.begin(), indices.end() );
  }

  {
  std::cout <<"Testing insertion of a range of segments\n";
  SDG sdg;
  std::vector< K::Segment_2 > segments;
  segments.push_back( K::Segment_2( points[0], points[1] ) );
  segments.push_back( K::Segment_2( points[1], points[2] ) );
  segments.push_back( K::Segment_2( points[2], points[3] ) );
  segments.push_back( K::Segment_2( points[3], points[4] ) );
  segments.push_back( K::Segment_2( points[4], points[0] ) );

  sdg.insert_segments( segments.begin(), segments.end() );
  }

  {
  std::cout <<"Testing insertion of a range of segments (pairs of points)\n";
  SDG sdg;
  std::vector< std::pair<K::Point_2, K::Point_2> > segments;
  segments.push_back( std::make_pair( points[0], points[1] ) );
  segments.push_back( std::make_pair( points[1], points[2] ) );
  segments.push_back( std::make_pair( points[2], points[3] ) );
  segments.push_back( std::make_pair( points[3], points[4] ) );
  segments.push_back( std::make_pair( points[4], points[0] ) );

  sdg.insert_segments( segments.begin(), segments.end() );
  }

  {
  std::cout <<"Testing insertion of a range of segments (Site_2)\n";
  SDG sdg;
  std::vector< Site_2 > segments;
  segments.push_back( Site_2::construct_site_2( points[0], points[1] ) );
  segments.push_back( Site_2::construct_site_2( points[1], points[2] ) );
  segments.push_back( Site_2::construct_site_2( points[2], points[3] ) );
  segments.push_back( Site_2::construct_site_2( points[3], points[4] ) );
  segments.push_back( Site_2::construct_site_2( points[4], points[0] ) );

  sdg.insert_segments( segments.begin(), segments.end() );
  }
}
