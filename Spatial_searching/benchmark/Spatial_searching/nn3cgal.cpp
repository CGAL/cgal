
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
typedef CGAL::Orthogonal_k_neighbor_search<TreeTraits_3,Distance> Neighbor_search_3;
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
  CGAL::set_binary_mode(data);
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
  int NN_number = 10;
  Points query_points_3, data_points_3;

  Point_3 p;
  Timer t;

  read(data_points_3, argv[1]);
  read(query_points_3, argv[2]);
 
  int bucketsize = (argc>3) ? boost::lexical_cast<int>(argv[3]) : 10;
  std::cerr << "bucketsize = "  << bucketsize <<std::endl;

  t.start();
  // Insert data points in the tree
  Tree_3 tree(data_points_3.begin(), data_points_3.end());
  tree.build();
  t.stop();
  std::cerr << "build " << t.time() << " sec\n";
  int items=0,leafs=0,internals=0;
  Points result(NN_number);

  t.reset();t.start();
  bool dump = true;
  for(Points::iterator it = query_points_3.begin(); it != query_points_3.end(); ++it) {
    // Initialize the search structure, and search all NN_number neighbors
    Neighbor_search_3 search(tree, *it, NN_number);
    int i=0;
    for (Neighbor_search_3::iterator it = search.begin(); it != search.end();it++, i++) {
      result[i] = it->first;
      if(dump){
        std::cerr << result[i].x()<<" "<<result[i].y()<<" "<<result[i].z() << std::endl;
		
      }
	}
	if(dump){
		tree.statistics(std::cerr);
		tree.print();
	}
    dump = false;
	//items += search.items_visited();
	//leafs += search.leaf_visited();
	//internals += search.internal_visited();
  }
  t.stop();
  std::cerr << "queries " << t.time() << " sec\n";
 // std::cerr << items <<" items "<<internals<< " internals "<<leafs<<  " leafs visited\n";
  std::cerr << "done\n";
  return 0;
}
