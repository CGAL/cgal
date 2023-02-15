#include<vector>

#include<boost/shared_ptr.hpp>

#include<CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include<CGAL/Polygon_2.h>
#include<CGAL/create_offset_polygons_2.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K ;

typedef K::Point_2                   Point ;
typedef CGAL::Polygon_2<K>           Polygon_2 ;
typedef CGAL::Straight_skeleton_2<K> Ss ;

typedef boost::shared_ptr<Polygon_2> PolygonPtr ;
typedef boost::shared_ptr<Ss> SsPtr ;

typedef std::vector<PolygonPtr> PolygonPtrVector ;

int main()
{
  std::ifstream in("data/A.poly");
  if (!in) return 1;
  int n;
  double x,y;
  Polygon_2 poly ;

  in >> n; // skip #polylines
  in >> n;
  for (int i=0; i<n; ++i)
  {
    in >> x >> y;
    poly.push_back( Point(x,y) ) ;
  }
  // we are only taking the outer contour (sufficient for this test)

  std::cout << poly.size() << "\n";

  SsPtr ss = CGAL::create_interior_straight_skeleton_2(poly);
  Ss ss_copy(*ss);

  double lOffset = 5 ;

  PolygonPtrVector offset_polygons = CGAL::create_offset_polygons_2<Polygon_2>(lOffset,ss_copy);

  assert(offset_polygons.size()==1);
  assert(offset_polygons.front()->size()>10);

  return 0;
}
