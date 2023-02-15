
#include <CGAL/Simple_cartesian.h>
#include <CGAL/K_neighbor_search.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Search_traits_3.h>
#include <CGAL/Timer.h>

#include <iostream>
#include <fstream>

#include <boost/lexical_cast.hpp>

typedef CGAL::Simple_cartesian<double> K;
typedef CGAL::Search_traits_3<K> TreeTraits_3;
typedef K::Point_3 Point_3;
typedef CGAL::Euclidean_distance<TreeTraits_3> Distance;
typedef CGAL::Sliding_midpoint<TreeTraits_3> Splitter;
typedef CGAL::Orthogonal_k_neighbor_search<TreeTraits_3,Distance,Splitter> Neighbor_search_3;
typedef Neighbor_search_3::Tree Tree_3;
typedef CGAL::Timer Timer;


typedef std::vector<Point_3> Points;


void read(Points& points, char* argv)
{
  int d, n;
  double x,y,z;
#if 0
  std::ifstream data(argv);
  data >> d >> n;
  assert(d == 3);
  points.reserve(n);
  while(data >> p){
    points.push_back(p);
  }
  data.close();

#else
 std::ifstream data(argv, std::ios::in | std::ios::binary);
  CGAL::IO::set_binary_mode(data);
  CGAL::read(data,d);
  CGAL::read(data,n);

  points.reserve(n);
  for(int i=0; i < n; i++){
    CGAL::read(data,x);
    CGAL::read(data,y);
    CGAL::read(data,z);
    points.push_back(Point_3(x,y,z));
  }
#endif
  data.close();
}

int main(int argc,char *argv[])
{
  Points query_points_3, data_points_3;

  Point_3 p;
  Timer t;

  read(data_points_3, argv[1]);
  read(query_points_3, argv[2]);

  int runs = (argc>3) ? boost::lexical_cast<int>(argv[3]) : 1;
  std::cerr << "runs = "  << runs <<std::endl;

  int bucketsize = (argc>4) ? boost::lexical_cast<int>(argv[4]) : 10;
  std::cerr << "bucketsize = "  << bucketsize <<std::endl;

  int NN_number = (argc>5) ? boost::lexical_cast<int>(argv[5]) : 10;
  std::cerr << "k = "  << NN_number <<std::endl;

  t.start();
  // Insert data points in the tree
  Tree_3 tree(data_points_3.begin(), data_points_3.end(),Splitter(bucketsize));
  tree.build();
  t.stop();
  std::cerr << "build " << t.time() << " sec\n";
  int items=0,leafs=0,internals=0;
  Points result(NN_number);
  tree.statistics(std::cerr);
  double sum=0;
  bool dump = true;


  for(int i = 0 ; i<runs; ++i){


    for(Points::iterator it = query_points_3.begin(); it != query_points_3.end(); ++it) {
      // Initialize the search structure, and search all NN_number neighbors
       t.reset();t.start();
      Neighbor_search_3 search(tree, *it, NN_number);
       t.stop();
      int i=0;
      for (Neighbor_search_3::iterator it = search.begin(); it != search.end();it++, i++) {
        result[i] = it->first;
        if(dump){
          std::cerr << result[i].x()<<" "<<result[i].y()<<" "<<result[i].z()<< std::endl;
        }
      }
      dump = false;
       sum += t.time();
       items+=search.items_visited();
       leafs+=search.leafs_visited();
       internals+=search.internals_visited();
    }

  }

  std::cerr << items <<" items\n";
  std::cerr << leafs <<" leaf\n";
  std::cerr << internals <<" internals visited\n";

  std::cerr<<std::endl << "total: " << sum << " sec\n";
  if(runs>1){
    std::cerr << "average: " << sum/runs << " sec\n";
  }
    std::cerr << "done\n";
  return 0;
}
