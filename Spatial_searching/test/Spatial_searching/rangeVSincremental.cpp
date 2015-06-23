#include <CGAL/Simple_cartesian.h>
#include <fstream>
#include <CGAL/Orthogonal_incremental_neighbor_search.h>
#include <CGAL/Search_traits_3.h>
#include <CGAL/Fuzzy_sphere.h>
#include <boost/function_output_iterator.hpp>

typedef CGAL::Simple_cartesian<double> K;
typedef CGAL::Search_traits_3<K> TreeTraits;
typedef CGAL::Orthogonal_incremental_neighbor_search<TreeTraits> NN_incremental_search;
typedef NN_incremental_search::iterator NN_iterator;
typedef NN_incremental_search::Tree Tree;
#include <vector>


typedef CGAL::Simple_cartesian<double>  K;
typedef K::Point_3 Point;

struct Inc {
  unsigned int * i;
  
  Inc(unsigned int& i)
    : i(&i)
  {}

  //  template <typename T>
  void operator()(const Point& t) const
  {
    std::cout << t << std::endl;
    ++(*i);
  }

};


int main(int argc, char* argv[]) {

  std::cout.precision(17);
  std::ifstream in(argv[1]);
  std::vector<Point> points;
  Point p;
  while(in >> p){
    points.push_back(p);
  }
  
  Tree tree;
  tree.insert(points.begin(), points.end());
 
  Point query(5,2,1);

  NN_incremental_search NN(tree, query);
  
  double sd = 0;
  std::cout <<  "The first 20 nearest neighbours with positive x-coord are: " << std::endl;
  NN_iterator it = NN.begin(); 
  for (int i=0; i<20; i++){
    sd = (*it).second;
    std::cout <<   (*it).first << "  at squared distance = " << (*it).second << std::endl;
    ++it;
  }

  CGAL::Fuzzy_sphere<TreeTraits> fs(query,sqrt(sd));
  std::vector<Point> result;
  tree.search(std::back_inserter(result), fs);
  for(int i=0; i < result.size(); i++){
    std::cout << result[i] << "   " << squared_distance(query, result[i]) << std::endl;
  }


  unsigned int count = 0;
  Inc inc(count);
  tree.search(boost::make_function_output_iterator(inc),fs);
  std::cout << count << std::endl;  
  

  return 0;
}
