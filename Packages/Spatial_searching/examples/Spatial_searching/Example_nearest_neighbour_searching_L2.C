
/*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
 * example_nearest_neighbour_searching_L2.C 
 * example using nearest_neighbour_iterator to illustrate
 * nearest neighbour searching using L2 distance
 *    
 *
 * Written by Hans Tangelder
 *            
 * 
\*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*/

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
#include <CGAL/Nearest_neighbour_L2.h>
#include <CGAL/Search_nearest_neighbour.h>

//#undef assert
//#define assert

#pragma hdrstop

#include <CGAL/point_generators_2.h>
#include <CGAL/algorithm.h>

int main() {

typedef CGAL::Cartesian<double> R;
typedef CGAL::Point_d<R> Point;
typedef Point::FT NT;
typedef CGAL::Plane_separator<NT> Separator;
typedef CGAL::Kd_tree_traits_point<Separator,Point> Traits;
typedef
CGAL::Nearest_neighbour_L2<Traits,CGAL::Search_nearest_neighbour>::iterator
NNN_Iterator;

{// start of scope of d

  int dim=2;
  int point_number=9801;
  int nearest_neighbour_number=5;
  int bucket_size=10;
  double eps=0.1;

  std::cout << " bucket_size=" << bucket_size << " eps=" << eps << std::endl;

  typedef std::list<Point> listd;
  listd lpt;

  std::cout << "Inserting evenly 9801 points  in the square (0,0)-(100,100) ...\n\n";
  for (int i=1; i<100; i++)
      for (int j=1; j<100; j++)
        {
          std::vector<double> vec(dim);
          vec[0]=double(i); vec[1]=double(j);
          Point p(dim, vec.begin(), vec.end());
          lpt.push_front(p);
        }


  typedef CGAL::Binary_search_tree<Traits> Tree;
  Tree d(lpt.begin(), lpt.end(), bucket_size);

  // test generating postscript file
  d.generate_postscript_file("test.ps",400,0,1);


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
        distance_to_q += (Q[i]-(*pli)[i])*(Q[i]-(*pli)[i]);
        distance_array[ii]=distance_to_q;
        ii++;
  };
  // t.stop();


  std::cout << "partial sort started" << std::endl;
  std::partial_sort(&distance_array[0],
                    &distance_array[2*nearest_neighbour_number],
                    &distance_array[point_number]);



  std::cout <<
  "comparison of approximate nearest neighbour distances and real distances"
  << std::endl;

  for (int i4=0; i4<nearest_neighbour_number; i4++) {
                 std::cout << "dist[" << i4 << "]=" << sqrt(result1[i4].second)
                 << std::endl;
                 std::cout << "distance_array[" << i4 << "]=" <<
                 sqrt(distance_array[i4]) << std::endl;
                 std::cout << "result[" << i4 << "]= " << *(result1[i4].first) << std::endl;
  };

  // compute the next k nearest neighbours

  std::cout << "testing iterator started computing next " << nearest_neighbour_number
  << " nearest neighbours" << std::endl;
  NNN_Iterator NNN_Iterator2=NNN_Iterator1;
  std::vector<Traits::Item_with_distance> result2(nearest_neighbour_number);

  for (int i3=0; i3 < nearest_neighbour_number; i3++) {
        result2[i3] = *NNN_Iterator2; ++NNN_Iterator2;
        std::cout << "dist[" << i3 << "]=" << sqrt(result2[i3].second)
        << std::endl;
        std::cout << "distance_array[" << nearest_neighbour_number+i3 << "]=" <<
        sqrt(distance_array[nearest_neighbour_number+i3]) << std::endl;
        std::cout << "result[" << i3 << "]= " << *(result2[i3].first) << std::endl;
  }

  std::cout << "testing iterator ready" << std::endl;

} // of scope of d

{ // start of scope of d2

  int dim=12;
  int point_number=5000;
  int nearest_neighbour_number=5;
  int bucket_size=10;
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



  typedef CGAL::Binary_search_tree<Traits> Tree;
  Tree d2(lpt.begin(), lpt.end(), bucket_size);
  
  // test generating postscript file
  // d2.generate_postscript_file("test.ps",400,0,1);


  std::vector<double> v(dim);
  for (int i2=0; i2<dim; i2++) v[i2]=1.0;
  v[0]=2.0;

  // query point
  Point Q(dim,v.begin(),v.end());
  std::cout << "query point is " << Q << std::endl;

  // illustrate use of container with copy

  std::cout << "started test copy" << std::endl;

  std::vector<Traits::Item_with_distance> result3(point_number);
  CGAL::Nearest_neighbour_L2<Traits,CGAL::Search_nearest_neighbour>
  NNN1(d2,Q,eps);

  std::copy(NNN1.begin(),NNN1.end(),result3.begin());

  std::cout << "copied approximate nearest neighbours are" << std::endl;

  for (int i3=0; i3 < nearest_neighbour_number; i3++) {
        std::cout << "result3[" << i3 << "]= " << *(result3[i3].first) << std::endl;
        std::cout << "with distance " << sqrt(result3[i3].second) << std::endl;
  }

  std::cout << "testing copy ready" << std::endl;

  // example of browsing
  std::cout << "started testing browsing with same input" << std::endl;

  CGAL::Nearest_neighbour_L2<Traits,CGAL::Search_nearest_neighbour> NNN2(d2,Q,eps);

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
}; // end of scope of d2

/*
double dummy;
std::cout << "Enter input to stop: \n" ;
std::cin >> dummy;
*/


return 0;
};



