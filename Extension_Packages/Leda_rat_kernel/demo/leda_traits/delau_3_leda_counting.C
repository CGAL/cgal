// this version counts kernel functor usages

// provide 3d kernel traits ...
#define CGAL_PROVIDE_LEDA_RAT_KERNEL_TRAITS_3
#define CGAL_NO_DEPRECATED_CODE

#include <CGAL/Cartesian.h>
#include <CGAL/Kernel_special.h>
#include <CGAL/kernel_counting_support.h>
#include <CEP/Leda_rat_kernel/leda_rat_kernel_traits.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <LEDA/integer.h>
#include <iostream>
#include <list>

#if defined(LEDA_NAMESPACE)
using namespace leda;
#endif

typedef CGAL::leda_rat_kernel_traits                    K1;
typedef CGAL::Do_counting<K1>                           Counting;
typedef CGAL::Kernel_special<K1, Counting >             K;

typedef CGAL::Delaunay_triangulation_3<K>               DT;
typedef K::Point_3                                      Point;


int main()
{
  DT T;

  // insertion of points
  int x,y,z;
  std::list<Point> L2;
  
  std::cout << "before construction !\n"; std::cout.flush();  
  
  for (z=0 ; z<10 ; z++)
    for (y=0 ; y<10 ; y++)
      for (x=0 ; x<10 ; x++) {
        leda_d3_rat_point pact(x,y,z,1);
        L2.push_back(pact); 
      }	  
  std::cout << "before insert !\n"; std::cout.flush();  	  
  T.insert(L2.begin(), L2.end());
  std::cout << "after insert !\n"; std::cout.flush();
  
  // check the result ...
  bool valid = T.is_valid(true);
  
  if (valid) std::cout << "valid result !\n";
  else std::cout << "NON - valid result !\n";
  
  std::cout << "Kernel functor usages:\n";
  std::cout << "----------------------\n";  
  Counting::print_counters(std::cout);
  
  return 0;
}
