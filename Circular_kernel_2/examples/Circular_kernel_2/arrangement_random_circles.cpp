#include <CGAL/Cartesian.h>
#include <CGAL/Random.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Algebraic_kernel_for_circles_2_2.h>
#include <CGAL/Circular_kernel_2.h>
#include <CGAL/Arr_circular_arc_traits.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Arr_naive_point_location.h>

typedef CGAL::Quotient<CGAL::MP_Float>                  NT;
typedef CGAL::Cartesian<NT>                             Linear_k;
typedef CGAL::Algebraic_kernel_for_circles_2_2<NT>      Algebraic_k;
typedef CGAL::Circular_kernel_2<Linear_k,Algebraic_k>   Circular_k;

typedef Circular_k::Point_2                             Point_2;
typedef Circular_k::Circle_2                            Circle_2;
typedef Circular_k::Circular_arc_2                      Circular_arc_2;
typedef std::vector<Circular_arc_2>                     ArcContainer;

typedef CGAL::Arr_circular_arc_traits<Circular_k>       Traits;

typedef CGAL::Arrangement_2<Traits>                     Arr;
typedef CGAL::Arr_naive_point_location<Arr>             Point_location;

int main(){

  CGAL::Random generatorOfgenerator;
  int random_seed = generatorOfgenerator.get_int(0, 123456);
  std::cout << "random_seed = " << random_seed << std::endl;
  CGAL::Random theRandom(random_seed);
  int random_max = 128;
  int random_min = -128;
  ArcContainer ac;
  int x, y;

  for (int i = 0; i < 10; i++) {
    x = theRandom.get_int(random_min,random_max);
    y = theRandom.get_int(random_min,random_max);
    ac.push_back( Circle_2( Point_2(x,y), x*x + y*y));
  }

  Arr  arr;
  Point_location _pl(arr);
  for (ArcContainer::const_iterator it=ac.begin(); it != ac.end(); ++it) {
    insert_curve(arr, *it, _pl);
  };

  return 0;
};
