// file: examples/Spatial_searching/General_neighbor_searching.C

#include <CGAL/Cartesian_d.h>
#include <CGAL/Manhattan_distance_iso_box_point.h>
#include <CGAL/K_neighbor_search.h>
#include <CGAL/Search_traits_2.h>

typedef CGAL::Cartesian_d<double> K;
typedef K::Point_d Point_d;
typedef K::Iso_box_d Iso_box_d;
typedef K TreeTraits;
typedef CGAL::Manhattan_distance_iso_box_point<TreeTraits> Distance;
typedef CGAL::K_neighbor_search<TreeTraits, Distance> Neighbor_search;
typedef Neighbor_search::Tree Tree;

int  main() {
  const int D = 2;
  const int N = 1;
  const int K = 1;

  std::list<Point_d> points;
  points.push_back(Point_d(0,0));

  Tree tree(points.begin(), points.end());

  Point_d pp(0,0);
  Point_d qq(1,1);
  Iso_box_d query(pp,qq);

  Distance tr_dist;
  Neighbor_search N1(tree, query, K, 10.0, false); // eps=10.0, nearest=false
  
  std::cout << "query = [0.1,0.2]^4 " << std::endl 
	    <<  K << " approximate furthest neighbors are: " << std::endl; 
  for (Neighbor_search::iterator it = N1.begin();it != N1.end();it++) { 
     std::cout << " d(q,fn)= " << tr_dist.inverse_of_transformed_distance(it->second) 
	       << " fn= " << it->first << std::endl; 
  } 
  return 0;
}
