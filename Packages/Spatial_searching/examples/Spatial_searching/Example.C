// example using nearest_neighbour_iterator for Linf

#include <CGAL/compiler_config.h>


#include <vector>
#include <numeric>
#include <cassert>

#include <iostream>
#include <CGAL/basic.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Point_d.h>
#include <CGAL/Binary_search_tree.h>
#include <CGAL/Kd_tree_traits_point.h>
#include <CGAL/Nearest_neighbour_Linf.h>
#include <CGAL/Search_nearest_neighbour.h>
#include <CGAL/Random.h>
#include <CGAL/Timer.h>

//#undef assert
//#define assert

#pragma hdrstop

// typedef CGAL::Point_d<R> Point;

#include <CGAL/point_generators_2.h>
#include <CGAL/algorithm.h>

int test_nearest_neighbour_Linf() {

typedef CGAL::Cartesian<double> R;
typedef CGAL::Point_d<R> Point;
typedef Point::FT NT;
typedef CGAL::Plane_separator<NT> Separator;
typedef CGAL::Kd_tree_traits_point<Separator,Point> Traits;
typedef CGAL::Nearest_neighbour_Linf<Traits,CGAL::Search_nearest_neighbour>::iterator NNN_Iterator;

{// start of scope of first test
  CGAL::Timer t;
  int dim=2;
  int point_number=1000;
  int nearest_neighbour_number=4;
  int bucket_size=10;
  double eps=0.1;

  std::cout << "test parameters: d=" << dim << " point_number=" << point_number << std::endl;
  std::cout << "nearest_neighbour_number=" << nearest_neighbour_number << " bucket_size="
  << bucket_size << " eps=" << eps << std::endl;

  typedef std::list<Point> listd;

  listd lpt;
  // add random points of dimension dim to lpt
  CGAL::Random Rnd;
  std::cout << "started testing nearest_neighbour_Linf" << std::endl;
  for (int i1=0; i1<point_number; i1++) {
        std::vector<double> vec(dim);
        for (int j=0; j<dim; j++) vec[j]=Rnd.get_double(-1.0,1.0);
        Point Random_point(dim, vec.begin(), vec.end());
        lpt.push_front(Random_point);
  }

  t.reset(); t.start();

  typedef CGAL::Binary_search_tree<Traits> Tree;
  Tree d(lpt.begin(), lpt.end(), bucket_size);
  t.stop();
  std::cout << "created binary search tree containing" << std::endl
  << point_number << " random points in the d-dim unit square in time " << t.time() <<
  " seconds " << std::endl;
  // d.statistics();

  d.generate_postscript_file("test.ps",400.0,0,1);

  // end of building binary search tree



  std::vector<double> v(dim);
  for (int i2=0; i2<dim; i2++) v[i2]=1.0;
  v[0]=2.0;

  // query point
  Point Q(dim,v.begin(),v.end());
  std::cout << "query point is " << Q << std::endl;

  NNN_Iterator NNN_Iterator1(d,Q,eps);

  std::vector<Traits::Item_with_distance> result1(nearest_neighbour_number);

  std::cout << "started testing copy_n" << std::endl;
  std::vector<Traits::Item_with_distance>::iterator it = result1.begin();
  CGAL::copy_n(NNN_Iterator1, nearest_neighbour_number, it);
  std::cout << "testing copy_n completed" << std::endl;

  // t.reset(); t.start();
  double distance_to_q;
  std::vector<double> distance_array(point_number);

  int ii=0;

  for (listd::iterator pli=lpt.begin(); pli != lpt.end(); pli++) {

        distance_to_q=0.0;

        for (int i=0; i<dim; i++)

        if ( fabs(Q[i]-(*pli)[i]) > distance_to_q )

        distance_to_q = fabs(Q[i]-(*pli)[i]);

        distance_array[ii]=distance_to_q;

        ii++;

  };
  // t.stop();

  /*
  std::cout << "computed distances in time " << t.time() <<

  " seconds " << std::endl;

  t.reset(); t.start();

  std::cout << "partial sort started" << std::endl; */

  std::partial_sort(&distance_array[0],

                    &distance_array[2*nearest_neighbour_number],

                    &distance_array[point_number]);

  /*

  t.stop();

  std::cout << "partial sort ready in time " << t.time() <<

  " seconds " << std::endl;

  */

  std::cout <<
  "comparison of approximate nearest neighbour distances and real distances"
  << std::endl;

  for (int i4=0; i4<nearest_neighbour_number; i4++) {
                 std::cout << "dist[" << i4 << "]=" << result1[i4].second
                 << std::endl;
                 std::cout << "distance_array[" << i4 << "]=" <<
                 distance_array[i4] << std::endl;
                 std::cout << "result[" << i4 << "]= " << *(result1[i4].first) << std::endl;
  };

  // compute the next k nearest neighbours

  std::cout << "number of items visited is " <<
  NNN_Iterator1.the_number_of_items_visited() << std::endl;

  std::cout << "testing iterator started computing next " << nearest_neighbour_number
  << " nearest neighbours" << std::endl;
  NNN_Iterator NNN_Iterator2=NNN_Iterator1;
  CGAL::swap<Traits,CGAL::Search_nearest_neighbour>(NNN_Iterator1, NNN_Iterator2);
  std::vector<Traits::Item_with_distance> result2(nearest_neighbour_number);

  for (int i3=0; i3 < nearest_neighbour_number; i3++) {
        result2[i3] = *NNN_Iterator2; ++NNN_Iterator2;
        std::cout << "dist[" << i3 << "]=" << result2[i3].second
        << std::endl;
        std::cout << "distance_array[" << nearest_neighbour_number+i3 << "]=" <<
        distance_array[nearest_neighbour_number+i3] << std::endl;
        std::cout << "result[" << i3 << "]= " << *(result2[i3].first) << std::endl;
  }

  std::cout << "number of items visited is " <<
  NNN_Iterator2.the_number_of_items_visited() << std::endl;

  std::cout << "testing iterator ready" << std::endl;

} // of scope of first test

{ // start of scope of second test
  CGAL::Timer t;
  int dim=12;
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
        std::vector<double> vec(dim);
        for (int j=0; j<dim; j++) vec[j]=Rnd.get_double(-1.0,1.0);
        Point Random_point(dim, vec.begin(), vec.end());
        lpt.push_front(Random_point);
  }

  t.reset(); t.start();

  typedef CGAL::Binary_search_tree<Traits> Tree;
  Tree d2(lpt.begin(), lpt.end(), bucket_size);
  t.stop();
  std::cout << "created binary search tree containing" << std::endl <<
  point_number << " random points in the d-dim unit square in time " << t.time() <<
  " seconds " << std::endl;
  // d2.statistics();

  // postscript stream not available yet in CGAL without using LEDA
  // d.generate_postscript_file("test.ps",12.5,0,1);

  // end of building binary search tree

  std::vector<double> v(dim);
  for (int i2=0; i2<dim; i2++) v[i2]=1.0;
  v[0]=2.0;

  // query point
  Point Q(dim,v.begin(),v.end());
  std::cout << "query point is " << Q << std::endl;

  // illustrate use of container with copy



  /*
  std::cout << "started test copy" << std::endl;

  std::vector<Traits::Item_with_distance> result3(point_number);
  CGAL::Nearest_neighbour_Linf<Traits> NNN1(d2,Q,eps);

  std::copy(NNN1.begin(),NNN1.end(),result3.begin());

  */

  /*
  if (NNN1.begin()==NNN1.end()) {
        std::cout << "side effect NNN.begin()==NNN.end() !!!" << std::endl;
  }
  else {
        std::cout << "expected side effect not detected ???" << std::endl;
  }*/


  /*
  std::cout << "copied approximate nearest neighbours are" << std::endl;

  for (int i3=0; i3 < nearest_neighbour_number; i3++) {
        std::cout << "result3[" << i3 << "]= " << *(result3[i3].first) << std::endl;
        std::cout << "with distance " << result3[i3].second << std::endl;
  }

  std::cout << "testing copy ready" << std::endl;
  */


  // example of browsing
  std::cout << "started testing browsing with same input" << std::endl;

  CGAL::Nearest_neighbour_Linf<Traits,CGAL::Search_nearest_neighbour> NNN2(d2,Q,eps);

  // define predicate class
  class GreaterThan0 {

        public:

                bool operator() (NNN_Iterator::value_type const result) const {

                return ( (*(result.first))[0] > 0.0);

        }

  };

  GreaterThan0 pred;

  NNN_Iterator first= NNN2.begin();
  NNN_Iterator last=NNN2.end();


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
}

int main() {
  test_nearest_neighbour_Linf();
  return 0;
}


