/*! \file
 * A Construction of a Planar_map and some point location queries.
 * A Construction of a Planar_map in the shape of the David star, some
 * point location queries of pre-coded points, and a point location
 * query of a point provided through the standard input.
 */

#include "short_names.h"

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Filtered_kernel.h>
#include <CGAL/Pm_segment_traits_2.h>
#include <CGAL/Pm_default_dcel.h>
#include <CGAL/Planar_map_2.h>
#include <iostream>
#include <iterator>
#include <algorithm>

typedef CGAL::Simple_cartesian<double>          Cartesian_kernel;
typedef CGAL::Filtered_kernel<Cartesian_kernel> Kernel;
typedef CGAL::Pm_segment_traits_2<Kernel>       Traits;
typedef Traits::Point_2                         Point_2;
typedef Traits::X_monotone_curve_2                       X_monotone_curve_2;
typedef CGAL::Pm_default_dcel<Traits>           Dcel;
typedef CGAL::Planar_map_2<Dcel,Traits>         Planar_map;
typedef Planar_map::Halfedge_handle             Halfedge_handle;
typedef Planar_map::Halfedge_const_handle       Halfedge_const_handle;
typedef Planar_map::Locate_type                 Locate_type;
typedef Planar_map::Ccb_halfedge_circulator     Ccb_halfedge_circulator;
typedef Planar_map::Ccb_halfedge_const_circulator
                                                Ccb_halfedge_const_circulator;

void print_point_locate(const Point_2 &, const Planar_map &);

int main()
{
  // Create an instance of a Planar_map:
  Planar_map pm;
  X_monotone_curve_2 cv[18];

  CGAL::set_ascii_mode(std::cout);

  Point_2 a1(6, 1), a2(1, 3), a3(4, 3), a4(8, 3), a5(11,3);
  Point_2 a6(3, 5), a7(9, 5), a8(1, 7), a9(4, 7), a10(8,7);
  Point_2 a11(11, 7), a12(6, 9);

  // Create the curves:
  cv[0]  = X_monotone_curve_2(a1, a3);
  cv[1]  = X_monotone_curve_2(a1, a4);
  cv[2]  = X_monotone_curve_2(a2, a3);
  cv[3]  = X_monotone_curve_2(a2, a6);
  cv[4]  = X_monotone_curve_2(a3, a6);
  cv[5]  = X_monotone_curve_2(a3, a4);
  cv[6]  = X_monotone_curve_2(a4, a5);
  cv[7]  = X_monotone_curve_2(a4, a7);
  cv[8]  = X_monotone_curve_2(a5, a7);
  cv[9]  = X_monotone_curve_2(a6, a8);
  cv[10] = X_monotone_curve_2(a6, a9);
  cv[11] = X_monotone_curve_2(a7, a10);
  cv[12] = X_monotone_curve_2(a7, a11);
  cv[13] = X_monotone_curve_2(a8, a9);
  cv[14] = X_monotone_curve_2(a9, a10);
  cv[15] = X_monotone_curve_2(a9, a12);
  cv[16] = X_monotone_curve_2(a10, a11);
  cv[17] = X_monotone_curve_2(a10, a12);

  std::cout << "The curves of the map :" << std::endl; 
  std::copy(&cv[0], &cv[18],
            std::ostream_iterator<X_monotone_curve_2>(std::cout, "\n"));
  std::cout << std::endl;

  // Insert the curves into the Planar_map:
  std::cout << "Inserting the curves to the map ... ";
  pm.insert(&cv[0], &cv[18]);
  std::cout << ((pm.is_valid()) ? "map valid!" : "map invalid!") << std::endl
            << std::endl;
  
  // Draw the map:
  std::cout << "    1  3 4  6  8 9 11		 " << std::endl;
  std::cout << "                                 "
               "                                " << std::endl;
  std::cout << "           a12            9             ";
  std::cout << "a1(" << a1 << ")  a2(" << a2 << ")  a3(" << a3 << ")";
  std::cout << std::endl;        
  std::cout << "           *  *           8             ";
  std::cout << "a4(" << a4 << ")  a5(" << a5 << ")  a6(" << a6 << ")";
  std::cout << std::endl;
  std::cout << "    a8===a9===a10==a11    7             ";
  std::cout << "a7(" << a7 << ")  a8(" << a8 << ")  a9(" << a9 << ")";
  std::cout << std::endl;
  std::cout << "     *  *       *  *      6             ";
  std::cout<< "a10(" << a10 << ")  a11(" << a11 << ")  a12(" << a12 << ")";
  std::cout << std::endl;
  std::cout << "      a6         a7       5" << std::endl;
  std::cout << "     *  *       *  *      4" << std::endl;
  std::cout << "    a2===a3===a4==a5      3" << std::endl;
  std::cout << "          *   *           2" << std::endl;
  std::cout << "           a1             1" << std::endl;

  // Find some point locations:
  print_point_locate(Point_2(6, 1), pm);
  print_point_locate(Point_2(2, 3), pm);
  print_point_locate(Point_2(6, 8), pm);
  std::cout << std::endl;

  std::cout << "Enter x coordinate for point location:" << std::endl;
  double x;
  std::cin >> x;
  std::cout << "Enter y coordinate for point location:" << std::endl;
  double y;
  std::cin >> y;
        
  print_point_locate(Point_2(x, y) ,pm);

  return 0;
}

/*! print_point_locate() locates the given point in the given Planar_map and
 * prints the computed location.
 * \param the point to be located
 * \param the Planar_map to operate on
 */
void print_point_locate(const Point_2 & p, const Planar_map & pm)
{
  Locate_type lt;
  Halfedge_const_handle edge = pm.locate(p, lt);
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

