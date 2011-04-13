#ifndef VOID_TEST
#include "numrep1.h"
#include <cassert>
#endif

#include <iostream>

#ifndef VOID_TEST
#include <CGAL/Pm_default_dcel.h>

#if TESTR == 1
#include <CGAL/Pm_segment_epsilon_traits.h>
#else
#include <CGAL/Pm_segment_exact_traits.h>
#endif

#if STRATEGY == 2
#include <CGAL/Pm_naive_point_location.h>
#else
#if STRATEGY == 3
#include <CGAL/Pm_walk_along_line_point_location.h>
#endif
#endif

#include <CGAL/Planar_map_2.h>

#include "numrep2.h"
typedef CGAL::Pm_default_dcel< Traits >     Dcel;
typedef CGAL::Planar_map_2< Dcel, Traits >  Planar_map;
typedef Traits::Point                       Point;
typedef Traits::X_curve                     Curve;
typedef Planar_map::Halfedge_handle         Halfedge_handle;
typedef Planar_map::Locate_type             Locate_type;
#endif // VOID_TEST

int main(int argc, char *argv[])
{
#ifndef VOID_TEST
#if STRATEGY == 3
  CGAL::Pm_walk_along_line_point_location<Planar_map> pl;
  Planar_map Pm(&pl);
#else
#if STRATEGY == 2
  CGAL::Pm_naive_point_location<Planar_map> pl;
  Planar_map Pm(&pl);
#else
  Planar_map Pm;
#endif
#endif
  
  int n; std::cin >> n;
  while (n--) {
    inputt x1, y1, x2, y2;
    std::cin >> x1 >> y1 >> x2 >> y2;
    Halfedge_handle hh = Pm.insert(Curve(Point(x1,y1),Point(x2,y2)));
    std::cout << "Inserted ("<< hh->curve() <<")"<<std::endl;
  }
  assert(Pm.is_valid());

  // General test
  // This used to be tst21, (Shai, Aug. 03, 2000)

  std::cout << "Faces" << Pm.number_of_faces() << std::endl;
  std::cout << "Halfedges" << Pm.number_of_halfedges() << std::endl;
  std::cout << "Vertices" << Pm.number_of_vertices() << std::endl;

  /*
  std::cout << "Clear planar map" << std::endl;
  Pm.clear();

  std::cout << "Faces" << Pm.number_of_faces() << std::endl;
  std::cout << "Halfedges" << Pm.number_of_halfedges() << std::endl;
  std::cout << "Vertices" << Pm.number_of_vertices() << std::endl;
  */

  // Locate test
  // This used to be tst22, (Shai, Aug. 03, 2000)
  std::cin >> n;
  while (n-->0) {
    inputt x,y;
    std::cin >> x >> y;
    Point p(x,y);
    std::cout << "Locate " << p <<std::endl;

    Locate_type lt;
     Halfedge_handle e = Pm.locate(p,lt);
    
    if (lt==Pm.UNBOUNDED_FACE) std::cout << "Unbounded face" << std::endl;
    else if (lt==Pm.FACE) std::cout << "Face that is left of " << e->curve() << std::endl;
    else if (lt==Pm.EDGE) std::cout << e->curve() << std::endl;
    //    else if (lt==Pm.VERTEX) std::cout << e.source().point() << std::endl;
    else if (lt==Pm.VERTEX) std::cout << e->target()->point() << std::endl; //changed 
    else std::cout << "Unknown locate type" << std::endl;
  }

  // Vertical Ray Shooting test
  // This used to be tst23, (Shai, Aug. 03, 2000)
  std::cin >> n;
  while (n-->0) {
    inputt x, y;
    std::cin >> x >> y;
    Point p(x,y);
    std::cout << "Vertical ray shoot " << p <<std::endl;

    Locate_type lt;
    Halfedge_handle e1 = Pm.vertical_ray_shoot(p,lt,true);
    
    std::cout << "Above" << std::endl;
    if (lt==Pm.UNBOUNDED_FACE) std::cout << "Unbounded face" << std::endl;
    //    else if (lt==Pm.EDGE) std::cout << e1->curve() << std::endl;
    else if (lt==Pm.EDGE) std::cout << e1->curve() << " (oriented toward " << e1->target()->point() << ")" << std::endl;
    //    else if (lt==Pm.VERTEX) std::cout << e1.source().point() << std::endl;
    else if (lt==Pm.VERTEX) std::cout << e1->target()->point() << std::endl; //changed
    else std::cout << "Unknown locate type" << std::endl;

    Halfedge_handle e2 = Pm.vertical_ray_shoot(p,lt,false);

    std::cout << "Below" << std::endl;
    if (lt==Pm.UNBOUNDED_FACE) std::cout << "Unbounded face" << std::endl;
    //    else if (lt==Pm.EDGE) std::cout << e2->curve() << std::endl;
    else if (lt==Pm.EDGE) std::cout << e2->curve() << " (oriented toward " << e2->target()->point() << ")" << std::endl;
    //    else if (lt==Pm.VERTEX) std::cout << e2.source().point() << std::endl;
    else if (lt==Pm.VERTEX) std::cout << e2->target()->point() << std::endl;
    else std::cout << "Unknown locate type" << std::endl;    
  }
#else
  usage();
#endif

  return 0;
}

