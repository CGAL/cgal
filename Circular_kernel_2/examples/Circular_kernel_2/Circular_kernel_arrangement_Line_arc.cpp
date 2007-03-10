#include <CGAL/Cartesian.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Gmpq.h>
#include <CGAL/Algebraic_kernel_for_circles_2_2.h>
#include <CGAL/intersections.h>
#include <CGAL/Circular_kernel_2.h>
#include <CGAL/Arr_line_arc_traits.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Arr_naive_point_location.h>
#include <CGAL/Random.h>

typedef CGAL::Quotient<CGAL::MP_Float>                      NT;
typedef CGAL::Cartesian<NT>                                 Linear_k;
typedef CGAL::Algebraic_kernel_for_circles_2_2<NT>          Algebraic_k;
typedef CGAL::Circular_kernel_2<Linear_k,Algebraic_k>       Circular_k;

//typedef Circular_k::Circular_arc_2                        Arc;
typedef Circular_k::Line_arc_2                              Arc_2;
typedef std::vector<Arc_2>                                  ArcContainer;

typedef CGAL::Arr_line_arc_traits<Circular_k>               Traits;

typedef CGAL::Arrangement_2<Traits>                         Arrangement;
typedef CGAL::Arr_naive_point_location<Arrangement>         Point_location;

typedef Circular_k::Point_2                           Point_2;
typedef Circular_k::Line_arc_2                        Line_arc_2;
int main(){

  CGAL::Random generatorOfgenerator;
  int random_seed = generatorOfgenerator.get_int(0, 123456);
  //std::cout << "random_seed = " << random_seed << std::endl;
  CGAL::Random theRandom(random_seed);
  int random_max = 5;
  int random_min = -5;
  for(int j = 0; j < 10 ; j++){

  ArcContainer ac;
  int x1;
  int y1;
  int x2;
  int y2;

    for(int i = 0; i < 50; i++){
      x1 = theRandom.get_int(random_min,random_max);
      y1 = theRandom.get_int(random_min,random_max);
      do{
   	x2 = theRandom.get_int(random_min,random_max);
   	y2 = theRandom.get_int(random_min,random_max);
      }while((x1 == x2) && (y1 == y2));


      std::cout << x1 << " "
   		<< y1 << " "
   		<< x2 << " "
   		<< y2 << std::endl;
      ac.push_back( Line_arc_2(Point_2(x1,y1), Point_2(x2,y2)));
    }

  //while(std::cin >> x1 >> y1 >> x2 >> y2){
   // ac.push_back( Line_arc_2(Point_2(x1,y1), Point_2(x2,y2)));
  //}

    {
      int i = 0;
      Arrangement arr;
      Point_location _pl(arr);
      for (ArcContainer::const_iterator it=ac.begin();
	   it != ac.end(); ++it) {
	std::cout << i++ << std::endl;
	//insert(arr,_pl,*it);
	insert_curve(arr, *it,_pl);
      }
    }
  }
  return 0;
};
