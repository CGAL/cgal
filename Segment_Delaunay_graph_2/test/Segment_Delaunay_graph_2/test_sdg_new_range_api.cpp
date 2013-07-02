#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Segment_Delaunay_graph_filtered_traits_2.h>
#include <CGAL/Segment_Delaunay_graph_2.h>
#include <CGAL/Gmpq.h>


typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Field_with_sqrt_tag MTag;
typedef CGAL::Integral_domain_without_division_tag EMTag;
typedef CGAL::Simple_cartesian<CGAL::Gmpq> EK;

typedef CGAL::Segment_Delaunay_graph_filtered_traits_without_intersections_2<K, MTag, EK, EMTag>  Gt;
typedef CGAL::Segment_Delaunay_graph_2<Gt>  SDG;
typedef SDG::Point_2 Point;
typedef SDG::Site_2 Site;


int main()
{
  std::vector<Point> points;
  points.push_back(Point(0,0));
  points.push_back(Point(0,1));
  points.push_back(Point(1,1));
  points.push_back(Point(2,3));
  points.push_back(Point(0,8));

  {
  SDG sdg;
  sdg.insert_points(points.begin(), points.end());
  }

  {
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
}
