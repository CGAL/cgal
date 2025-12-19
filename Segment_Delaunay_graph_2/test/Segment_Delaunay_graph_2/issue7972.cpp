
#include <CGAL/Segment_Delaunay_graph_2.h>
#include <CGAL/Segment_Delaunay_graph_traits_2.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

typedef CGAL::Exact_predicates_exact_constructions_kernel K;
typedef CGAL::Segment_Delaunay_graph_traits_2<K>          Gt;
typedef CGAL::Segment_Delaunay_graph_2<Gt>                SDG2;

int main()
{
  std::vector<CGAL::Segment_2<K>> segments;
  segments.emplace_back(CGAL::Point_2<K>(0.0, 0.0),CGAL::Point_2<K>(1.0, 0.0));

  SDG2 delaunay;
  delaunay.insert_segments(segments.begin(), segments.end());

  return 0;
}
