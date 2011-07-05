
#include <CGAL/Simple_cartesian.h>
#include <cassert>
#include <CGAL/point_generators_2.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/Random.h>
#include <CGAL/Search_traits_2.h>
#include <CGAL/Search_traits_3.h>
#include <CGAL/Orthogonal_incremental_neighbor_search.h>
#include <CGAL/algorithm.h>



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

  typedef std::vector<Point>               Vector;



template <class S>
struct Splitter_test {


  typedef CGAL::Orthogonal_incremental_neighbor_search<TreeTraits,Distance,S> Orthogonal_incremental_neighbor_search;
  typedef typename Orthogonal_incremental_neighbor_search::iterator NN_iterator;
  typedef typename Orthogonal_incremental_neighbor_search::Tree Tree;
  typedef typename Orthogonal_incremental_neighbor_search::Point_with_transformed_distance Point_with_transformed_distance;
  
  struct less_second {
    bool operator()(const typename Orthogonal_incremental_neighbor_search::Point_with_transformed_distance& p1,
		    const typename Orthogonal_incremental_neighbor_search::Point_with_transformed_distance& p2) const
    {
      return p1.second < p2.second;
    }
  };
  
  /*
    std::ostream& operator<<(std::ostream& os, const Orthogonal_incremental_neighbor_search::Point_with_transformed_distance& pair)
    {
    os << pair.first << " | " << pair.second;
    return os;
    }
  */

  bool operator()(Vector points, Point query) {

    Vector points2;
    
    Tree t(points.begin(), points.end());
    
    std::vector<typename Orthogonal_incremental_neighbor_search::Point_with_transformed_distance> result; 
    
    Orthogonal_incremental_neighbor_search oins(t, query, 0.0, true);
    
    typename Orthogonal_incremental_neighbor_search::iterator it = oins.begin();
    Point_with_transformed_distance pd = *it;
    points2.push_back(pd.first);
    if(CGAL::squared_distance(query,pd.first) != pd.second){
      std::cout << CGAL::squared_distance(query,pd.first) << " != " << pd.second << std::endl;
    }
    assert(CGAL::squared_distance(query,pd.first) == pd.second);
    it++;
    for(; it != oins.end();it++){
      Point_with_transformed_distance qd = *it;
      assert(pd.second <= qd.second);
      pd = qd;
      points2.push_back(pd.first);
      if(CGAL::squared_distance(query,pd.first) != pd.second){
	std::cout << CGAL::squared_distance(query,pd.first) << " != " << pd.second << std::endl;
      }
      assert(CGAL::squared_distance(query,pd.first) == pd.second);
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
  CGAL::copy_n( g, 1000, std::back_inserter(points));
    
  g++;
  Point query = *g;
  {
    std::cout << "Testing Sliding_midpoint" << std::endl;
    Splitter_test<CGAL::Sliding_midpoint<TreeTraits> > st;
    st(points,query);
  }
  {
    std::cout << "Testing Sliding_fair" << std::endl;
    Splitter_test<CGAL::Sliding_fair<TreeTraits> > st;
    st(points,query);
  }
  
  {
    std::cout << "Testing Fair" << std::endl;
    Splitter_test<CGAL::Fair<TreeTraits> > st;
    st(points,query);
  }
  
  {
    std::cout << "Testing Median_of_rectangle" << std::endl;
    Splitter_test<CGAL::Median_of_rectangle<TreeTraits> > st;
    st(points,query);
  }

  
  {
    std::cout << "Testing Midpoint_of_rectangle" << std::endl;
    Splitter_test<CGAL::Midpoint_of_rectangle<TreeTraits> > st;
    st(points,query);
  }
  
  {
    std::cout << "Testing Midpoint_of_max_spread" << std::endl;
    Splitter_test<CGAL::Midpoint_of_max_spread<TreeTraits> > st;
    st(points,query);
  }
  {
    std::cout << "Testing Median_of_max_spread" << std::endl;
    Splitter_test<CGAL::Median_of_max_spread<TreeTraits> > st;
    st(points,query);
  }

  std::cout << "done" << std::endl;
  return 0;
}  

