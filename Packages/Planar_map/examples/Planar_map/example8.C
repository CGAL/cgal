//example8

#include <iostream>

#define CGAL_NO_PM_DEFAULT_POINT_LOCATION
#define NAIVE_POINT_LOCATION
#define CARTESIAN	//  Cartesian or Homogeneous
#define ROBUSTNESS

#if defined(CGAL_USE_LEDA) && defined(ROBUSTNESS)
#include <CGAL/leda_rational.h>
#endif

#ifdef CARTESIAN
#include <CGAL/Cartesian.h>
#else
#include <CGAL/Homogeneous.h>
#endif

#include <CGAL/Pm_straight_exact_traits.h>

#include <CGAL/Pm_default_dcel.h>
#include <CGAL/Planar_map_2.h>

#include <CGAL/IO/Straight_2_stream.h>
#include <CGAL/IO/Planar_map_iostream.h>

#ifdef NAIVE_POINT_LOCATION
#include <CGAL/Pm_naive_point_location.h>
#endif

using namespace CGAL;

#if defined(CGAL_USE_LEDA) && defined(ROBUSTNESS)
#ifdef CARTESIAN
typedef Cartesian<leda_rational>	Coord_t;
#else
typedef Homogeneous<leda_rational>	Coord_t;
#endif
#else //defined(CGAL_USE_LEDA) && defined(ROBUSTNESS)
#ifdef CARTESIAN
typedef Cartesian<double>		Coord_t;
#else
typedef Homogeneous<double>		Coord_t;
#endif
#endif

typedef Pm_straight_exact_traits<Coord_t>	Traits;
typedef Traits::Bounding_box			Bounding_box;
typedef Traits::Point				Point;
typedef Traits::X_curve				X_curve;
typedef Traits::X_bounded_curve			Segment;
typedef Pm_default_dcel<Traits>			Dcel;
typedef Planar_map_2<Dcel,Traits>		PM;
typedef PM::Halfedge_handle			Halfedge_handle;

int main()
{
	
  Bounding_box bbox(Point(-10,-10),Point(10,10));
  Traits traits(bbox);

  // creating an instance of Planar_map_Bbox_2<Dcel,Traits>
#ifndef NAIVE_POINT_LOCATION   
  PM pm(traits);
#else	
  Pm_naive_point_location< PM > naive_pl;
  PM pm(traits, &naive_pl, NULL);
#endif
  
  X_curve cv[6];
  int i;
  
  set_ascii_mode(std::cout);
  Point a1(1, 1), a2(1, 0), a3(0, 0), a4(0, 1), a5(1,4,2) ;
  
  /*
	a5 
        /\  
       /  \
   a4  ----  a1
      |    | 
      |    | 
      |    |
   a3  ----  a2
 */        
	
	
  // those curves are about to enter to pm
  cv[0] = X_curve(Segment(a1, a2));
  cv[1] = X_curve(Segment(a2, a3));
  cv[2] = X_curve(Segment(a3, a4));
  cv[3] = X_curve(Segment(a4, a5));
  cv[4] = X_curve(Segment(a5, a1));
  cv[5] = X_curve(Segment(a1, a4));
  
  std::cout << "the curves of the map :" << std::endl; 
  for (i = 0; i < 6; i++)
    std::cout << cv[i] << std::endl;
  
  std::cout << std::endl;
  
  
  // insert the five curves to the map
  std::cout << "inserting the curves to the map..." << std::endl;
  Halfedge_handle e[6];  
  
  e[0]=pm.insert(cv[0]);
  
  for (i = 1; i < 4; i++)
    {
      e[i]=pm.insert_from_vertex(cv[i],e[i-1]->target(), true);
    }
  
  e[4]=pm.insert_at_vertices(cv[4],e[0]->source(),e[3]->target() );
  
  e[5]=pm.insert_at_vertices(cv[5],e[0]->source(),e[2]->target() );
	
  /* 
         O            
     e3 /\  e4
       /  \
      O----O
      | e5 | 
   e2 |    | e0
      |    |
      O----O
	e1 
 */
	
  // check the validity of the map
  std::cout << "check map validity... ";
  if (pm.is_valid())
    std::cout << "map valid!" << std::endl;
  else
    std::cout << "map invalid!" << std::endl;
  std::cout << std::endl << "Printing map: " << std::endl;
  
  std::cout << pm ;
  
  return 0;  
}
