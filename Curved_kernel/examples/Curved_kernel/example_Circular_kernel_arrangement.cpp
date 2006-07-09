//#include <fstream>

#include <CGAL/basic.h>
#include <CGAL/Cartesian.h>
#include <CGAL/point_generators_2.h>

#include <CGAL/MP_Float.h>
#include <CGAL/Gmpq.h>

#include <CGAL/Algebraic_kernel_2_2.h>

#include <CGAL/intersections.h>

#include <CGAL/Circular_kernel.h>
#include <CGAL/Arr_circular_arc_traits.h>
#include <CGAL/Arr_circular_arc_traits_tracer.h>

#include <CGAL/Arrangement_2.h>
#include <CGAL/Arr_naive_point_location.h>



#include <CGAL/Random.h>


// typedef CGAL::MP_Float                                          NT;
typedef CGAL::Gmpq                                          NT;
// typedef CGAL::Lazy_exact_nt<CGAL::Gmpq>                     NT;
// typedef CGAL::Lazy_exact_nt<CGAL::Quotient<CGAL::MP_Float> >    NT;
// typedef CGAL::Lazy_exact_nt<CGAL::Quotient<CGAL::Gmpz> >    NT;
typedef CGAL::Cartesian<NT>                                 Linear_k;

typedef CGAL::Algebraic_kernel_for_circles_2_2<NT>          Algebraic_k;
typedef CGAL::Circular_kernel_2<Linear_k,Algebraic_k>       Circular_k;
// typedef CGAL::Circular_kernel<Linear_k>                        Circular_k;

typedef Circular_k::Circular_arc_2                            Arc;
typedef std::vector<Arc>                                    ArcContainer;

#ifndef CGAL_CURVED_KERNEL_DEBUG
typedef CGAL::Arr_circular_arc_traits<Circular_k>                  Traits;
#else
typedef CGAL::Arr_circular_arc_traits<Circular_k>                  Traits0;
typedef CGAL::Circular_arc_traits_tracer<Traits0>            Traits;
#endif

//typedef Traits::Point_2                             Point_2;
typedef Traits::Curve_2                             Conic_arc_2;
typedef CGAL::Arrangement_2<Traits>                 Pmwx;
typedef CGAL::Arr_naive_point_location<Pmwx>        Point_location;

typedef Traits::X_monotone_curve_2                          X_monotone_curve_2;
typedef Circular_k::Point_2                    Point_2;
typedef Circular_k::Circle_2                    Circle_2;
typedef Circular_k::Circular_arc_2              Circular_arc_2;
int main(){
  
  CGAL::Random generatorOfgenerator;
  int random_seed = generatorOfgenerator.get_int(0, 123456);
  std::cout << "random_seed = " << random_seed << std::endl;
  CGAL::Random theRandom(random_seed);
  int random_max = 128;
  int random_min = -128;
  ArcContainer ac;
  int x;
  int y;
  for(int i = 0; i < 10; i++){
    x = theRandom.get_int(random_min,random_max);
    y = theRandom.get_int(random_min,random_max);
    ac.push_back( Circle_2( Point_2(x,y), x*x + y*y));
  }

  Pmwx _pm;
  Point_location _pl(_pm);
  for (ArcContainer::const_iterator it=ac.begin();
       it != ac.end(); ++it) {
    //insert(_pm,_pl,*it);
    insert_curve(_pm,*it,_pl);
      };
  
  return 0;
};

