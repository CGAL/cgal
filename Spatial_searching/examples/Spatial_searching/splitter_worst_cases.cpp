#include <CGAL/Simple_cartesian.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Search_traits_2.h>

typedef CGAL::Simple_cartesian<double>                Kernel;
typedef CGAL::Search_traits_2<Kernel>                 Traits_2;
typedef Kernel::Point_2                               Point_2;
typedef CGAL::Sliding_midpoint<Traits_2>              Sliding_midpoint;
typedef CGAL::Median_of_rectangle<Traits_2>           Median_of_rectangle;
typedef CGAL::Euclidean_distance<Traits_2>            Distance;
typedef CGAL::Orthogonal_k_neighbor_search<Traits_2,Distance,Sliding_midpoint>     Neighbor_search_sliding;
typedef CGAL::Orthogonal_k_neighbor_search<Traits_2,Distance,Median_of_rectangle>  Neighbor_search_median;
typedef Neighbor_search_sliding::Tree                 Tree_sliding;
typedef Neighbor_search_median::Tree                  Tree_median;

typedef std::vector<Point_2>                          Points;

int main()
{
  Points sliding_worst_case;
  for (int i = 0 ,j = 1; i < 10 ; ++i , j *= 2){
    sliding_worst_case.push_back(Point_2(((double)i)/10 , 0));
    sliding_worst_case.push_back(Point_2( (double)j , 0));    
  }

  Sliding_midpoint sliding(10);
  Median_of_rectangle median(10);

  Tree_sliding tree1(sliding_worst_case.begin(), sliding_worst_case.end() , sliding);
  tree1.build();

  std::cout << "Worst case tree for Sliding midpoint and Midpoint of max spread : "<<std::endl;
  tree1.statistics(std::cout);
  tree1.clear();
  std::cout<<std::endl<<"Same data with median splitter:"<<std::endl;

  Tree_median tree2(sliding_worst_case.begin(), sliding_worst_case.end() , median );
  tree2.statistics(std::cout);
  tree2.clear();

  Points median_worst_case;
  for(int i = 0 ; i < 19 ; ++i){
    median_worst_case.push_back(Point_2( 0 , i));
  }
  median_worst_case.push_back(Point_2(20,0));

 
  Tree_median tree3(median_worst_case.begin() , median_worst_case.end() , median);

  tree3.build();
  std::cout <<std::endl<< "Worst case tree for Median of rectangle, Median of max spread : "<<std::endl;
  tree3.statistics(std::cout);
  tree3.clear();

  std::cout<<std::endl<<"Same data with midpoint splitter:"<<std::endl;

  Tree_sliding tree4(median_worst_case.begin() , median_worst_case.end() , sliding);

  tree4.build();
  tree4.statistics(std::cout);
  tree4.clear();
  return 0;
}