
#include <CGAL/Simple_cartesian.h>
#include <cassert>
#include <CGAL/point_generators_2.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/point_generators_d.h>
#include <CGAL/Random.h>
#include <CGAL/Search_traits_2.h>
#include <CGAL/Search_traits_3.h>
#include <CGAL/Orthogonal_incremental_neighbor_search.h>
#include <CGAL/Incremental_neighbor_search.h>
#include <CGAL/algorithm.h>



#ifdef TWO
typedef CGAL::Simple_cartesian<double>         K;
typedef K::Point_2                       Point;
typedef CGAL::Creator_uniform_2<double,Point>  Creator;
typedef CGAL::Random_points_in_disc_2<Point,Creator> Random_points;
typedef CGAL::Search_traits_2<K> TreeTraits;
#else
typedef CGAL::Simple_cartesian<double>         K;
typedef K::Point_3                       Point;
typedef CGAL::Creator_uniform_3<double,Point>  Creator;
typedef CGAL::Random_points_in_sphere_3<Point,Creator> Random_points;
typedef CGAL::Search_traits_3<K> TreeTraits;
#endif



typedef std::vector<Point>               Vector;
//typedef CGAL::Orthogonal_incremental_neighbor_search<TreeTraits> Orthogonal_incremental_neighbor_search;
typedef CGAL::Incremental_neighbor_search<TreeTraits> Orthogonal_incremental_neighbor_search;
typedef Orthogonal_incremental_neighbor_search::iterator NN_iterator;
typedef Orthogonal_incremental_neighbor_search::Tree Tree;
typedef Orthogonal_incremental_neighbor_search::Point_with_transformed_distance Point_with_transformed_distance;

struct less_second {
  bool operator()(const Orthogonal_incremental_neighbor_search::Point_with_transformed_distance& p1,
		  const Orthogonal_incremental_neighbor_search::Point_with_transformed_distance& p2) const
  {
    return p1.second < p2.second;
  }
};


std::ostream& operator<<(std::ostream& os, const Orthogonal_incremental_neighbor_search::Point_with_transformed_distance& pair)
{
  os << pair.first << " | " << pair.second;
  return os;
}

int 
main() {
  const int N = 1000;
  Vector points;
  Vector points2;
  Random_points g( 150.0);

  // We enforce IEEE double precision as we compare a distance 
  // in a register with a distance in memory
  CGAL::Set_ieee_double_precision pfr;

  CGAL::copy_n( g, 1000, std::back_inserter(points));

  Tree t(points.begin(), points.end());
  g++;
  Point query = *g;

  std::vector<Orthogonal_incremental_neighbor_search::Point_with_transformed_distance> result(N); 

  Orthogonal_incremental_neighbor_search oins(t, query, 0.0, true);

  Orthogonal_incremental_neighbor_search::iterator it = oins.begin();
  Point_with_transformed_distance pd = *it;
  points2.push_back(pd.first);
  if(CGAL::squared_distance(query,pd.first) != pd.second){
    std::cout << "different distances: " << CGAL::squared_distance(query,pd.first) << " != " << pd.second << std::endl;
  }

  assert(CGAL::squared_distance(query,pd.first) == pd.second);
  it++;
  for(; it != oins.end();it++){
    Point_with_transformed_distance qd = *it;
    assert(pd.second <= qd.second);
    pd = qd;
    points2.push_back(pd.first);
    if(CGAL::squared_distance(query,pd.first) != pd.second){
      std::cout  << "different distances: " << CGAL::squared_distance(query,pd.first) << " != " << pd.second << std::endl;
    }
    assert(CGAL::squared_distance(query,pd.first) == pd.second);
  }


  std::sort(points.begin(),points.end());
  std::sort(points2.begin(),points2.end());
  assert(points.size() == points2.size());
  for(unsigned int i = 0; i < points.size(); i++){
    if(points[i] != points2[i]){
      std::cout << i << "  " << points[i] << "   " << points2[i] << std::endl;
    }
  }
  assert(points == points2);
  

  std::cout << "done" << std::endl;

  return 0;
}
  


