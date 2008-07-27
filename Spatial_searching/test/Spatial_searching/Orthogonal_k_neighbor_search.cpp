// file          : test/Spatial_searching/Orthogonal_k_neighbor_search.C

#include <CGAL/Cartesian.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/Search_traits_2.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <set>

typedef CGAL::Cartesian<double> K;
typedef K::Point_2 Point;
typedef CGAL::Random_points_in_square_2<Point> Random_points_iterator;
typedef CGAL::Counting_iterator<Random_points_iterator> N_Random_points_iterator;
typedef CGAL::Search_traits_2<K> TreeTraits;
typedef CGAL::Orthogonal_k_neighbor_search<TreeTraits> Neighbor_search;
typedef Neighbor_search::Tree Tree;


void 
search(bool nearest)
{
  const unsigned int N = 1000;
  const double cube_side_length = 1.0;
  // generator for random data points in the square ( (-1,-1), (1,1) ) 
  Random_points_iterator rpit( cube_side_length);

  std::vector<Point> points(N_Random_points_iterator(rpit,0),
			    N_Random_points_iterator(N));
 
  
  Tree tree(points.begin(), points.end());
  Point query(0,0);
  
  Neighbor_search search(tree, query, N/2 , 0.0, nearest);
 
  std::vector<Point> result, diff; 
   // report the N/2 furthest neighbors and their distance

  //std::copy(search.begin(), search.end(), std::back_inserter(result));
  for(Neighbor_search::iterator nit = search.begin();
      nit != search.end();
      nit++){
    result.push_back(nit->first);
  }

  std::sort(points.begin(), points.end());
  std::sort(result.begin(), result.end());
  std::set_difference(points.begin(), points.end(),
		      result.begin(), result.end(),
		      std::back_inserter(diff));

  std::cout << "|result| = " << result.size() << "  |diff| = " << diff.size() << std::endl;
  double sep_dist = (nearest)?0:3 * cube_side_length;
  {
    for(std::vector<Point>::iterator it = result.begin();
	it != result.end();
	it++){
      double dist = CGAL::squared_distance(query, *it);
      if(nearest){
	if(dist > sep_dist) sep_dist = dist;
      } else {
	if(dist < sep_dist) sep_dist = dist;
      }
    }
  }
  // the other points must be further/closer than min_dist
  {
    for(std::vector<Point>::iterator it = diff.begin();
	it != diff.end();
	it++){
      double dist = CGAL::squared_distance(query, *it);
      if(nearest){
	if(dist < sep_dist){
	  std::cout << "Error: Point " << *it << " at distance " << dist << "  <  " << sep_dist << std::endl;
	}
      } else {
	if(dist > sep_dist){
	  std::cout << "Error: Point " << *it << " at distance " << dist << "  >  " << sep_dist  << std::endl;
	}
      }
    }
  }
}

int main()
{
  search(true);
  search(false);
  std::cout << "done" << std::endl;
  return 0;
}
