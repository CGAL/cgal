// example using nearest_neighbour_iterator for Linf
// illustrating approximate k-nearest and 
// approximate next nearest neighbour search


#include <CGAL/compiler_config.h>

#include <CGAL/basic.h>

#include <vector>
#include <list>
#include <numeric>
#include <cassert>

#include <iostream>

#include <CGAL/Cartesian_d.h>
#include <CGAL/Binary_search_tree.h>
#include <CGAL/Kd_tree_traits_point.h>
#include <CGAL/Nearest_neighbour_Linf.h>
#include <CGAL/Search_nearest_neighbour.h>
#include <CGAL/Random.h>
#include <CGAL/Timer.h>




#include <CGAL/point_generators_2.h>
#include <CGAL/algorithm.h>

int test_nearest_neighbour_Linf() {

typedef CGAL::Cartesian_d<double> R;
typedef CGAL::Point_d<R> Point;

typedef CGAL::Kernel_traits<Point>::Kernel K;
typedef K::FT NT;
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
        double v[] = {Rnd.get_double(-1.0,1.0), Rnd.get_double(-1.0,1.0)};
        Point Random_point(dim,v,v+2);
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

  //  d.generate_postscript_file("test.ps",400.0,0,1);

  // end of building binary search tree


  /*
  std::vector<double> v(dim);
  */
  
  // Query point;
  double v[] = {2.0,1.0};
  Point Q(2,v,v+2);
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


