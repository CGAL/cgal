// file: test/Spatial_searching/K_neighbor_search_manhattan_distance_isobox_point.C

#include <CGAL/Cartesian.h>
#include <cassert>
#include <CGAL/Manhattan_distance_iso_box_point.h>
#include <CGAL/K_neighbor_search.h>
#include <CGAL/Search_traits_2.h>
#include <CGAL/Search_traits_adapter.h>
#include "Point_with_info.h"

typedef CGAL::Cartesian<double> Kernel;
typedef Kernel::Point_2 Point_d;
typedef Kernel::Iso_rectangle_2 Iso_box_d;
typedef CGAL::Search_traits_2<Kernel> TreeTraits;
typedef CGAL::Manhattan_distance_iso_box_point<TreeTraits> Distance;
typedef CGAL::K_neighbor_search<TreeTraits, Distance> Neighbor_search;
//typdefs for Point_with_info
typedef Point_with_info_helper<Point_d>::type                                           Point_with_info;
typedef Point_property_map<Point_d>                                                     Ppmap;
typedef CGAL::Search_traits_adapter<Point_with_info,Ppmap,TreeTraits>                   Traits_with_info;
typedef CGAL::Distance_adapter <Point_with_info,Ppmap,Distance>                         Distance_adapter;
typedef CGAL::K_neighbor_search<Traits_with_info, Distance_adapter>                   Neighbor_search_with_info;

const unsigned int K = 8;


template <class K_search>
void run(std::list<Point_d> points)
{
  typename K_search::Tree tree(
    boost::make_transform_iterator(points.begin(),Create_point_with_info<typename K_search::Point_d>()),
    boost::make_transform_iterator(points.end(),Create_point_with_info<typename K_search::Point_d>())
  );

  Point_d pp(2,3);
  Point_d qq(4,5);
  Iso_box_d query(pp,qq);

  Distance tr_dist;
  K_search N1(tree, query, K); // eps=10.0, nearest=false
  
  for (typename K_search::iterator it = N1.begin();it != N1.end();it++) { 
    assert( get_point(it->first) == points.front());
    points.pop_front();
  } 
  std::cout << "done" << std::endl;  
}

int  main() {
  

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

  run<Neighbor_search>(points);
  run<Neighbor_search_with_info>(points);

  return 0;
}
