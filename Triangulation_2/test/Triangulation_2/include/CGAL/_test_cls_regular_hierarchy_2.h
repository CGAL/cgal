//#include <CGAL/point_generators_2.h>
//#include <CGAL/iterator.h>
#include <CGAL/Random.h>
#include <CGAL/algorithm.h>
#include <CGAL/_test_cls_regular_triangulation_2.h>

template <class Rh>
void
_test_cls_regular_hierarchy_2( const Rh & )
{
  typedef Rh  Regular_hierarchy;
  typedef typename Regular_hierarchy::Weighted_point   Weighted_point;
  typedef typename Regular_hierarchy::Bare_point       Bare_point;

  _test_cls_regular_triangulation_2( Regular_hierarchy());

  //int nn = 500;
  int nn=100;
  std::cout << " insertion of " << nn << "  points" << std::endl;
  Regular_hierarchy rh;

  CGAL::Random rand;
  // std::ofstream output("data"); CGAL::IO::set_binary_mode(output);
  for(int i = 0; i < nn ; i++){
    Bare_point p( rand.get_double(), rand.get_double());
    Weighted_point wp(p, (rand.get_double())*100);
    //output <<  wp ;
    //std::cerr << i << " " << std::endl ;
    rh.insert(wp);
    rh.is_valid();
  }
  std::cerr << std::endl;
  rh.is_valid(true);

//    std::ifstream input("data"); CGAL::IO::set_binary_mode(input);
//    Weighted_point wp;
//    int inr = 0;
//    while(input) {
//      inr++;
//      input >> wp;
//      std::cerr << inr << " wpoint " << wp.point() << " " <<  wp.weight()
//                <<std::endl;
//      rh.insert(wp);
//    }
//    rh.is_valid(true);


//    int sign = 1;
//    int nn = 50;
//    for(int i = 0; i <nn; i++) {
//      for(int j = 0; j<nn; j++) {
//        sign= -sign;
//        Bare_point p( i,j);
//        Weighted_point wp(p,500 + sign*500);
//        rh.insert(wp);
//      }
//    }

//rh.is_valid(true);

  std::cout << "  location" << std::endl;
  rh.locate(Weighted_point(0.,0.));

  std::cout <<  "  removal of all points" << std::endl;
  while( rh.number_of_vertices() > 0) {
    rh.remove(rh.finite_vertices_begin());
  }
  return;
}




