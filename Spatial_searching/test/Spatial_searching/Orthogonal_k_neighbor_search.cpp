// file          : test/Spatial_searching/Orthogonal_k_neighbor_search.C

#include <CGAL/Simple_cartesian.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/Search_traits_2.h>
#include <CGAL/Search_traits_adapter.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include "Point_with_info.h"
#include <set>

typedef CGAL::Simple_cartesian<double> K;
typedef K::Point_2 Point;
typedef CGAL::Random_points_in_square_2<Point>                          Random_points_iterator;
typedef CGAL::Counting_iterator<Random_points_iterator>                 N_Random_points_iterator;
typedef CGAL::Search_traits_2<K>                                        TreeTraits;
typedef CGAL::Orthogonal_k_neighbor_search<TreeTraits>                  Neighbor_search;
//typdefs fo Point_with_info
typedef Point_with_info_helper<Point>::type                                                             Point_with_info;
typedef Point_property_map<Point>                                                  Ppmap;
typedef CGAL::Search_traits_adapter<Point_with_info,Ppmap,TreeTraits>                        Traits_with_info;
typedef CGAL::Distance_adapter <Point_with_info,Ppmap,Neighbor_search::Distance>   Distance_adapter;
typedef CGAL::Orthogonal_k_neighbor_search<Traits_with_info,Distance_adapter>                         Neighbor_search_with_info;

template <class K_search>
bool search(bool nearest)
{
  const unsigned int N = 1000;
  const double cube_side_length = 1.0;
  // generator for random data points in the square ( (-1,-1), (1,1) )
  Random_points_iterator rpit( cube_side_length);

  std::vector<Point> points(N_Random_points_iterator(rpit,0),
                            N_Random_points_iterator(N));


  typename K_search::Tree tree(
    boost::make_transform_iterator(points.begin(),Create_point_with_info<typename K_search::Point_d>()),
    boost::make_transform_iterator(points.end(),Create_point_with_info<typename K_search::Point_d>())
  );
  Point query(0,0);

  K_search search(tree, query, N/2 , 0.0, nearest);

  std::vector<Point> result, diff;
  // report the N/2 furthest neighbors and their distance

  //std::copy(search.begin(), search.end(), std::back_inserter(result));
  for(typename K_search::iterator nit = search.begin();
      nit != search.end();
      nit++){
    result.push_back(get_point(nit->first));
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
  bool res=true;
  // the other points must be further/closer than min_dist
  {
    for(std::vector<Point>::iterator it = diff.begin();
        it != diff.end();
        it++){
      double dist = CGAL::squared_distance(query, *it);
      if(nearest){
        if(dist < sep_dist){
          std::cout << "Error: Point " << *it << " at distance " << dist << "  <  " << sep_dist << std::endl;
          res=false;
        }
      } else {
        if(dist > sep_dist){
          std::cout << "Error: Point " << *it << " at distance " << dist << "  >  " << sep_dist  << std::endl;
          res=false;
        }
      }
    }
  }
  return res;
}

int main()
{
  bool res=true;
  res&=search<Neighbor_search>(true);
  res&=search<Neighbor_search>(false);
  res&=search<Neighbor_search_with_info>(true);
  res&=search<Neighbor_search_with_info>(false);
  std::cout << "done" << std::endl;
  return res ? 0 : 1;
}
