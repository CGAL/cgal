#ifndef VOID_TEST
#include "numrep1.h"
#include <cassert>
#endif //VOID_TEST

#include <iostream>

#ifndef VOID_TEST
#include <CGAL/Pm_default_dcel.h>
#include <CGAL/Planar_map_2.h>
#include <CGAL/IO/Straight_2_stream.h>
#include <CGAL/Pm_straight_exact_traits.h>

#if STRATEGY == 2
#include <CGAL/Pm_naive_point_location.h>
#else
#if STRATEGY == 3
#include <CGAL/Pm_walk_along_line_point_location.h>
#endif
#endif

#include "numrep2.h"
#include <CGAL/Pm_dynamic_open_bounding_box.h>

typedef CGAL::Pm_default_dcel< Traits > Dcel;
typedef CGAL::Planar_map_2< Dcel, Traits>  Planar_map;
typedef Traits::Point Point;
typedef Traits::X_curve Curve;
typedef Traits::X_bounded_curve Segment;
typedef Traits::X_target_unbounded_curve Ray;
typedef Traits::X_unbounded_curve Line;
typedef Planar_map::Halfedge_handle Halfedge_handle;
typedef Planar_map::Locate_type Locate_type;
#endif //VOID_TEST

int main(int argc, char *argv[])
{
#ifndef VOID_TEST
#if TESTR != 3
  Traits tr;
#else
  Traits tr(2);// work with integers, using GMP NT.
#endif
  CGAL::Pm_dynamic_open_bounding_box<Planar_map> bb; 
#if STRATEGY == 3
  CGAL::Pm_walk_along_line_point_location<Planar_map> pl;
  Planar_map Pm(tr,&pl,&bb);
#else
#if STRATEGY == 2
  CGAL::Pm_naive_point_location<Planar_map> pl;
  Planar_map Pm(tr,&pl,&bb);
#else
  CGAL::Pm_default_point_location<Planar_map> pl;
  Planar_map Pm(tr,&pl,&bb);
#endif
#endif
  
  int n; std::cin >> n;
  while (n--) {
    inputt x1,y1,x2,y2;
    std::cin >> x1 >> y1 >> x2 >> y2;
    Halfedge_handle e = Pm.insert(Segment(Point(x1,y1),Point(x2,y2)));
    std::cout << "\nInserted Segment("<< e->curve() <<")"<<std::flush;
  }
  std::cin >> n;
  while (n--) {
    inputt x1, y1, x2, y2;
    std::cin >> x1 >> y1 >> x2 >> y2;
    Halfedge_handle e = Pm.insert(Ray(Point(x1,y1),Point(x2,y2)));
    std::cout << "\nInserted Ray("<< e->curve() <<")"<<std::flush;
  }
  std::cin >> n;
  while (n--) {
    inputt a,b,c;
    std::cin >> a >> b >> c;
    Halfedge_handle e = Pm.insert(Line(a,b,c));
    std::cout << "\nInserted Line("<< e->curve() <<")"<<std::flush;
  }
  // point location queries
  std::cin >> n;
  while (n--) {
    inputt a,b;
    std::cin >> a >> b;
    Locate_type lt;
    Halfedge_handle e = Pm.locate(Point(a,b),lt);
    std::cout << "\nLocate("<< a << "," << b << ")=";
    if (lt==Planar_map::VERTEX) 
      std::cout << "VERTEX " << e->target()->point() << std::endl;
    else if (lt==Planar_map::EDGE)
      std::cout << "EDGE " << e->curve()
	//		<< " oriented toward " << e->target()->point() 
		<< std::endl;
    else if (lt==Planar_map::UNBOUNDED_FACE)
      std::cout << "UNBOUNDED_FACE" << std::endl;
    else if (lt==Planar_map::FACE)
      std::cout << "FACE (to left of) " << e->curve()
	//		<< " oriented toward " << e->target()->point() 
		<< std::endl;
    else if (lt==Planar_map::UNBOUNDED_VERTEX) 
      std::cout << "UNBOUNDED_VERTEX " << e->target()->point() << std::endl;
    else if (lt==Planar_map::UNBOUNDED_EDGE)
      std::cout << "UNBOUNDED_EDGE " << e->curve()
	//		<< " oriented toward " << e->target()->point() 
		<< std::endl;
    else
      CGAL_assertion(
                     lt==Planar_map::VERTEX||lt==Planar_map::EDGE||
                     lt==Planar_map::FACE || lt==Planar_map::UNBOUNDED_FACE ||
                     lt==Planar_map::UNBOUNDED_VERTEX||
                     lt==Planar_map::UNBOUNDED_EDGE);
  }
#else //VOID_TEST
  usage();
#endif //VOID_TEST

  return 0;
}






