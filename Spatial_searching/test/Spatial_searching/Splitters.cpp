
#include <CGAL/Simple_cartesian.h>
#include <cassert>
#include <CGAL/point_generators_2.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/Random.h>
#include <CGAL/Search_traits_2.h>
#include <CGAL/Search_traits_3.h>
#include <CGAL/Search_traits_adapter.h>
#include <CGAL/Orthogonal_incremental_neighbor_search.h>
#include <CGAL/algorithm.h>
#include "Point_with_info.h"



#ifdef TWO
typedef CGAL::Simple_cartesian<double>         K;
typedef K::Point_2                       Point;
typedef CGAL::Creator_uniform_2<double,Point>  Creator;
typedef CGAL::Random_points_in_disc_2<Point,Creator> Random_points;
typedef CGAL::Search_traits_2<K> TreeTraits;
typedef CGAL::Euclidean_distance<TreeTraits> Distance;
#else
typedef CGAL::Simple_cartesian<double>         K;
typedef K::Point_3                       Point;
typedef CGAL::Creator_uniform_3<double,Point>  Creator;
typedef CGAL::Random_points_in_sphere_3<Point,Creator> Random_points;
typedef CGAL::Search_traits_3<K> TreeTraits;
typedef CGAL::Euclidean_distance<TreeTraits> Distance;
#endif
typedef Point_with_info_helper<Point>::type                                             Point_with_info;
typedef Point_property_map<Point>                                                       Ppmap;
typedef CGAL::Search_traits_adapter<Point_with_info,Ppmap,TreeTraits>                   Traits_with_info;
typedef CGAL::Distance_adapter <Point_with_info,Ppmap,Distance>                         Distance_adapter;

  typedef std::vector<Point>               Vector;



template <class Traits,class S,class Distance>
struct Splitter_test {


  typedef CGAL::Orthogonal_incremental_neighbor_search<Traits,Distance,S> Orthogonal_incremental_neighbor_search;
  typedef typename Orthogonal_incremental_neighbor_search::iterator NN_iterator;
  typedef typename Orthogonal_incremental_neighbor_search::Tree Tree;
  typedef typename Orthogonal_incremental_neighbor_search::Point_with_transformed_distance Point_with_transformed_distance;

  bool operator()(Vector points, Point query) {

    Vector points2;

    Tree t(
      boost::make_transform_iterator(points.begin(),Create_point_with_info<typename Traits::Point_d>()),
      boost::make_transform_iterator(points.end(),Create_point_with_info<typename Traits::Point_d>())
    );

    Orthogonal_incremental_neighbor_search oins(t, query, 0.0, true);

    typename Orthogonal_incremental_neighbor_search::iterator it = oins.begin();
    Point_with_transformed_distance pd = *it;
    points2.push_back(get_point(pd.first));
    if(CGAL::squared_distance(query,get_point(pd.first)) != pd.second){
      std::cout << CGAL::squared_distance(query,get_point(pd.first)) << " != " << pd.second << std::endl;
    }
    assert(CGAL_IA_FORCE_TO_DOUBLE(CGAL::squared_distance(query,get_point(pd.first))) == pd.second);
    it++;
    for(; it != oins.end();it++){
      Point_with_transformed_distance qd = *it;
      assert(pd.second <= qd.second);
      pd = qd;
      points2.push_back(get_point(pd.first));
      if(CGAL_IA_FORCE_TO_DOUBLE(CGAL::squared_distance(query,get_point(pd.first))) != pd.second){
        std::cout << CGAL::squared_distance(query,get_point(pd.first)) << " != " << pd.second << std::endl;
      }
      assert(CGAL_IA_FORCE_TO_DOUBLE(CGAL::squared_distance(query,get_point(pd.first))) == pd.second);
    }
    std::sort(points.begin(),points.end());
    std::sort(points2.begin(),points2.end());
    assert(points.size() == points2.size());
    assert(points == points2);
    return true;
  }
};

int
main() {

  // We enforce IEEE double precision as we compare a distance
  // in a register with a distance in memory
  CGAL::Set_ieee_double_precision pfr;

  Vector points;
  Random_points g( 150.0);
  std::copy_n( g, 1000, std::back_inserter(points));

  g++;
  Point query = *g;
  {
    std::cout << "Testing Sliding_midpoint" << std::endl;
    Splitter_test<TreeTraits,CGAL::Sliding_midpoint<TreeTraits>,Distance > st;
    st(points,query);
    Splitter_test<Traits_with_info,CGAL::Sliding_midpoint<Traits_with_info>,Distance_adapter > sti;
    sti(points,query);
  }
  {
    std::cout << "Testing Sliding_fair" << std::endl;
    Splitter_test<TreeTraits,CGAL::Sliding_fair<TreeTraits>,Distance > st;
    st(points,query);
    Splitter_test<Traits_with_info,CGAL::Sliding_fair<Traits_with_info>,Distance_adapter > sti;
    sti(points,query);
  }

  {
    std::cout << "Testing Fair" << std::endl;
    Splitter_test<TreeTraits,CGAL::Fair<TreeTraits>,Distance > st;
    st(points,query);
    Splitter_test<Traits_with_info,CGAL::Fair<Traits_with_info>,Distance_adapter > sti;
    sti(points,query);
  }

  {
    std::cout << "Testing Median_of_rectangle" << std::endl;
    Splitter_test<TreeTraits,CGAL::Median_of_rectangle<TreeTraits>,Distance > st;
    st(points,query);
    Splitter_test<Traits_with_info,CGAL::Median_of_rectangle<Traits_with_info>,Distance_adapter > sti;
    sti(points,query);
  }


  {
    std::cout << "Testing Midpoint_of_rectangle" << std::endl;
    Splitter_test<TreeTraits,CGAL::Midpoint_of_rectangle<TreeTraits>,Distance > st;
    st(points,query);
    Splitter_test<Traits_with_info,CGAL::Midpoint_of_rectangle<Traits_with_info>,Distance_adapter > sti;
    sti(points,query);
  }

  {
    std::cout << "Testing Midpoint_of_max_spread" << std::endl;
    Splitter_test<TreeTraits,CGAL::Midpoint_of_max_spread<TreeTraits>,Distance > st;
    st(points,query);
    Splitter_test<Traits_with_info,CGAL::Midpoint_of_max_spread<Traits_with_info>,Distance_adapter > sti;
    sti(points,query);
  }
  {
    std::cout << "Testing Median_of_max_spread" << std::endl;
    Splitter_test<TreeTraits,CGAL::Median_of_max_spread<TreeTraits>,Distance > st;
    st(points,query);
    Splitter_test<Traits_with_info,CGAL::Median_of_max_spread<Traits_with_info>,Distance_adapter > sti;
    sti(points,query);
  }

  std::cout << "done" << std::endl;
  return 0;
}

