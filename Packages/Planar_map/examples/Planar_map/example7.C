// file: examples/Planar_map/example7.C

/*! \file
 * A construction of a house-like-shaped Planar_map.
 * In this example a house-like-shaped Planar_map that uses the three point
 * location strategies listed below is constructed.
 * 1. Trapezoidal Decomposition
 * 2. Naive
 * 3. Walk
 */

#include "short_names.h"

#include <CGAL/Cartesian.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Quotient.h>
#include <CGAL/Pm_segment_traits_2.h>
#include <CGAL/Pm_default_dcel.h>
#include <CGAL/Planar_map_2.h>

#include <CGAL/Pm_trapezoid_ric_point_location.h>
#include <CGAL/Pm_walk_along_line_point_location.h>
#include <CGAL/Pm_naive_point_location.h>

#include <iterator>
#include <algorithm>

typedef CGAL::Quotient<CGAL::MP_Float>                  Number_type;
typedef CGAL::Cartesian<Number_type>                    Kernel;
typedef CGAL::Pm_segment_traits_2<Kernel>               Traits;
typedef CGAL::Pm_default_dcel<Traits>                   Dcel;
typedef CGAL::Planar_map_2<Dcel,Traits>                 Planar_map;
typedef CGAL::Pm_trapezoid_ric_point_location<Planar_map>     Trap_point_location;
typedef CGAL::Pm_naive_point_location<Planar_map>       Naive_point_location;
typedef CGAL::Pm_walk_along_line_point_location<Planar_map>
                                                        Walk_point_location;
typedef Planar_map::Locate_type                         Locate_type;
typedef Traits::Point_2                                 Point_2;
typedef Traits::X_monotone_curve_2                      X_monotone_curve_2;
typedef Planar_map::Halfedge_const_handle               Halfedge_const_handle;
typedef Planar_map::Ccb_halfedge_const_circulator
                                                Ccb_halfedge_const_circulator;

void print_point_locate(const Point_2 & p, const Planar_map * pm)
{
  Locate_type lt;
  Halfedge_const_handle edge = pm->locate(p, lt);
  Ccb_halfedge_const_circulator curr, first;
  std::cout << "The location of point " << p << " is of type ";
  switch (lt) {
    case Planar_map::VERTEX :
      std::cout << "VERTEX" << std::endl
                << "The vertex is: (" << edge->target()->point() << ")" 
                << std::endl;
      break;

    case Planar_map::UNBOUNDED_VERTEX :
      std::cout << "UNBOUNDED_VERTEX" << std::endl
                << "The vertex is: (" << edge->target()->point() << ")" 
                << std::endl;
      break;

    case Planar_map::EDGE :
      std::cout << "EDGE" << std::endl
                << "The edge is: {(" << edge->source()->point()
                << ")->(" << edge->target()->point() << ")}" << std::endl;
      break;
      
    case Planar_map::UNBOUNDED_EDGE :
      std::cout << "UNBOUNDED_EDGE" << std::endl
                << "The edge is: {(" << edge->source()->point()
                << ")->(" << edge->target()->point() << ")}" << std::endl;
      break;

    case Planar_map::FACE :
      first = Ccb_halfedge_const_circulator(edge);
      curr = first;
      std::cout << "FACE" << std::endl
                << "The face is: [" << "(" << curr->target()->point() << ")";
      for (++curr; curr != first; ++curr)
        std::cout << ", (" << curr->target()->point() << ")";
      std::cout << "]" << std::endl;
      break;

    case Planar_map::UNBOUNDED_FACE :
      std::cout << "UNBOUNDED_FACE" << std::endl;
      break;
  }
}

int main(int argc, char * argv[])
{
  Point_2 a1(1, 1), a2(1, 0), a3(0, 0), a4(0, 1), a5(1, 4, 2);

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

  // Create the curves:
  X_monotone_curve_2 cv[6];
  cv[0] = X_monotone_curve_2(a1, a2);
  cv[1] = X_monotone_curve_2(a2, a3);
  cv[2] = X_monotone_curve_2(a3, a4);
  cv[3] = X_monotone_curve_2(a4, a5);
  cv[4] = X_monotone_curve_2(a5, a1);
  cv[5] = X_monotone_curve_2(a1, a4);

  CGAL::set_ascii_mode(std::cout);
  std::cout << "The curves of the map :" << std::endl; 
  std::copy(&cv[0], &cv[6],
            std::ostream_iterator<X_monotone_curve_2>(std::cout, "\n"));
  std::cout << std::endl;

#define PL_TRAP         0
#define PL_NAIVE        1
#define PL_WALK         2
#define NUM_PLS         3
  Planar_map * pms[NUM_PLS];
  
  // Trap:
  Trap_point_location trap_pl;
  Planar_map pm_trap(&trap_pl);
  pms[PL_TRAP] = &pm_trap;

    // Naive:
  Naive_point_location naive_pl;  
  Planar_map pm_naive(&naive_pl);
  pms[PL_NAIVE] = &pm_naive;

  // Walk:
  Walk_point_location walk_pl;
  Planar_map pm_walk(&walk_pl);
  pms[PL_WALK] = &pm_walk;

  for (int i = 0; i < NUM_PLS; i++) {
    pms[i]->insert(&cv[0], &cv[6]);
    std::cout << ((pms[i]->is_valid()) ? "map valid!" : "map invalid!")
              << std::endl;
    print_point_locate(Point_2(1, 4, 2), pms[i]);
    print_point_locate(Point_2(1, 2, 2), pms[i]);
    print_point_locate(Point_2(1, 5, 4), pms[i]);
    print_point_locate(Point_2(1, 0, 2), pms[i]);
    print_point_locate(Point_2(2, 0), pms[i]);
    std::cout << std::endl;
  }
  
  return 0;  
}
