// example using nearest_neighbour_iterator for Linf
// testing approximate browsing of all items in a search tree
// testing browsing using a predicate to find the first 
// nearest neighbour with a positive value in the first dimension


#include <CGAL/compiler_config.h>

#include <CGAL/basic.h>

#include <vector>
#include <list>
#include <numeric>
#include <cassert>

#include <iostream>

#include <CGAL/Cartesian_d.h>
//#include <CGAL/Point_d.h>
#include <CGAL/Binary_search_tree.h>
#include <CGAL/Kd_tree_traits_point.h>
#include <CGAL/Nearest_neighbour_Linf.h>
#include <CGAL/Search_nearest_neighbour.h>
#include <CGAL/Random.h>
#include <CGAL/Timer.h>

//#undef assert
//#define assert

#pragma hdrstop



#include <CGAL/point_generators_2.h>
#include <CGAL/algorithm.h>

typedef CGAL::Cartesian_d<double> R;
typedef CGAL::Point_d<R> Point;

typedef CGAL::Kernel_traits<Point>::Kernel K;
typedef K::FT NT;
typedef CGAL::Plane_separator<NT> Separator;
typedef CGAL::Kd_tree_traits_point<Separator,Point> MyTraits;
typedef CGAL::Nearest_neighbour_Linf<MyTraits,CGAL::Search_nearest_neighbour>::iterator NNN_Iterator;

int test_nearest_neighbour_Linf() {


{ // start of scope of second test
  CGAL::Timer t;
  const int dim=12;
  int point_number=5;
  int nearest_neighbour_number=5;
  int bucket_size=1;
  double eps=0.1;

  std::cout << "test parameters: d=" << dim << " point_number=" << point_number << std::endl;
  std::cout << "nearest_neighbour_number=" << nearest_neighbour_number << " bucket_size="
  << bucket_size << " eps=" << eps << std::endl;

  typedef std::list<Point> listd;

  listd lpt;
  // add random points of dimension dim to lpt
  CGAL::Random Rnd;
  // std::cout << "started tstrandom()" << std::endl;
  for (int i1=0; i1<point_number; i1++) {
	    double v[dim];
		for (int i2=0; i2<dim; i2++) v[i2]=Rnd.get_double(-1.0,1.0);
        Point Random_point(dim,v,v+dim);
        lpt.push_front(Random_point);
  }

  t.reset(); t.start();

  typedef CGAL::Binary_search_tree<MyTraits> Tree;
  Tree d2(lpt.begin(), lpt.end(), bucket_size);
  t.stop();
  std::cout << "created binary search tree containing" << std::endl <<
  point_number << " random points in the d-dim unit square in time " << t.time() <<
  " seconds " << std::endl;
  d2.statistics();

  
  // d2.generate_postscript_file("test2.ps",400.0,0,1);

  // end of building binary search tree

  double vq[dim];
  for (int i2=0; i2<dim; i2++) vq[i2]=1.0;
  vq[0]=2.0;

  // query point
  Point Q(dim,vq,vq+dim);
  std::cout << "query point is " << Q << std::endl;

  // illustrate use of container with copy



  
  std::cout << "started test copy" << std::endl;

  std::vector<MyTraits::Item_with_distance> result3(point_number);
  CGAL::Nearest_neighbour_Linf<MyTraits,CGAL::Search_nearest_neighbour> NNN1(d2,Q,eps);

  NNN_Iterator begin, end;
  begin = NNN1.begin();
  end = NNN1.end();
  std::copy(begin,end,result3.begin());

  

  
  if (NNN1.begin()==NNN1.end()) {
        std::cout << "side effect NNN.begin()==NNN.end() !!!" << std::endl;
  }
  else {
        std::cout << "expected side effect not detected ???" << std::endl;
  }


  
  std::cout << "copied approximate nearest neighbours are" << std::endl;

  for (int i3=0; i3 < nearest_neighbour_number; i3++) {
        std::cout << "result3[" << i3 << "]= " << *(result3[i3].first) << std::endl;
        std::cout << "with distance " << result3[i3].second << std::endl;
  }

  std::cout << "testing copy ready" << std::endl;
  


  // example of browsing
  std::cout << "started testing browsing with same input" << std::endl;

  CGAL::Nearest_neighbour_Linf<MyTraits,CGAL::Search_nearest_neighbour> NNN2(d2,Q,eps);



  // define predicate class
  class GreaterThan0 {

        public:

                bool operator() (NNN_Iterator::value_type const result) const {

                return ( (*(result.first))[0] > 0.0);

        }

  };

  GreaterThan0 pred;

  NNN_Iterator first;
  first = NNN2.begin();
  NNN_Iterator last;
  last=NNN2.end();


  // find_if did not work, but
  // a copy of code of find_if  does
  while (first != last && !pred(*first)) ++first;

  if (last != first)  {
        std::cout << "first positive neighbour is " << (*(*first).first)
        << std::endl;
  }
  else  {
        std::cout << "no positive neighbour found" << std::endl;
  };
  std::cout << "testing browsing ready" << std::endl;
}; // end of scope of second test

return 0;
};

int main() {
  test_nearest_neighbour_Linf();
  
  /*
  double dummy;
  std::cout << "Enter input to stop: \n" ;
  std::cin >> dummy;
  */

  return 0;
};


