// file: test/Spatial_searching/K_neighbor_search_manhattan_distance_isobox_point.C

#include <CGAL/Cartesian.h>
#include <cassert>
#include <CGAL/Manhattan_distance_iso_box_point.h>
#include <CGAL/K_neighbor_search.h>
#include <CGAL/Search_traits_2.h>

typedef CGAL::Cartesian<double> K;
typedef K::Point_2 Point_d;
typedef K::Iso_rectangle_2 Iso_box_d;
typedef CGAL::Search_traits_2<K> TreeTraits;
typedef CGAL::Manhattan_distance_iso_box_point<TreeTraits> Distance;
typedef CGAL::K_neighbor_search<TreeTraits, Distance> Neighbor_search;
typedef Neighbor_search::Tree Tree;

int  main() {
  const unsigned int K = 8;

  std::list<Point_d> points;
  points.push_back(Point_d(3,4));
  points.push_back(Point_d(2,2));
  points.push_back(Point_d(3,1));
  points.push_back(Point_d(6,2));
  points.push_back(Point_d(8,4));
  points.push_back(Point_d(8,6));
  points.push_back(Point_d(9,6));
  points.push_back(Point_d(8,8));
  points.push_back(Point_d(9,1));
  points.push_back(Point_d(10,1));
  points.push_back(Point_d(10,9));


  Tree tree(points.begin(), points.end());

  Point_d pp(2,3);
  Point_d qq(4,5);
  Iso_box_d query(pp,qq);

  Distance tr_dist;
  Neighbor_search N1(tree, query, K); // eps=10.0, nearest=false
  
  for (Neighbor_search::iterator it = N1.begin();it != N1.end();it++) { 
    assert( it->first == points.front());
    points.pop_front();
  } 
  std::cout << "done" << std::endl;
  return 0;
}
