/*! \file
 * A construction of a house-like-shaped Planar_map.
 * In this example a house-like-shaped Planar_map that uses the naive point
 * location strategy is constructed utilising various edge-insertion
 * operations, that is, insert(), insert_from_vertex(), and
 * insert_at_vertices(). The naive point location strategy suits better
 * these functions. The planar map is printed out to the standard output using
 * the provided streaming feature.
 */

#include "short_names.h"

#include <CGAL/Cartesian.h>
#include <CGAL/Quotient.h>
#include <CGAL/Pm_segment_traits_2.h>
#include <CGAL/Pm_default_dcel.h>
#include <CGAL/Planar_map_2.h>
#include <CGAL/IO/Pm_file_writer.h>
#include <CGAL/IO/write_pm.h>
#include <CGAL/Pm_naive_point_location.h>
#include <iostream>
#include <iterator>
#include <algorithm>

typedef CGAL::Quotient<float>                   Number_type;
typedef CGAL::Cartesian<Number_type>            Kernel;
typedef CGAL::Pm_segment_traits_2<Kernel>       Traits;
typedef Traits::Point_2                         Point_2;
typedef Traits::X_monotone_curve_2              X_monotone_curve_2;
typedef CGAL::Pm_default_dcel<Traits>           Dcel;
typedef CGAL::Planar_map_2<Dcel,Traits>         Planar_map;
typedef CGAL::Pm_file_writer<Planar_map>        Pm_writer;

int main()
{
  // Create an instance of a planar map:
  CGAL::Pm_naive_point_location<Planar_map> naive_pl;
  Planar_map pm(&naive_pl);
  Pm_writer verbose_writer(std::cout, pm, true);
  X_monotone_curve_2 cv[6];
  int i;

  CGAL::set_ascii_mode(std::cout);

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
  cv[0] = X_monotone_curve_2(a1, a2);
  cv[1] = X_monotone_curve_2(a2, a3);
  cv[2] = X_monotone_curve_2(a3, a4);
  cv[3] = X_monotone_curve_2(a4, a5);
  cv[4] = X_monotone_curve_2(a5, a1);
  cv[5] = X_monotone_curve_2(a1, a4);

  std::cout << "The curves of the map :" << std::endl; 
  std::copy(&cv[0], &cv[6],
            std::ostream_iterator<X_monotone_curve_2>(std::cout, "\n"));
  std::cout << std::endl;

  // Insert the curves into the Planar_map:
  std::cout << "Inserting the curves into the map ... ";
  Planar_map::Halfedge_handle e[6];
  e[0] = pm.insert(cv[0]);
  for (i = 1; i < 4; i++)
    e[i] = pm.insert_from_vertex(cv[i], e[i-1]->target());
  e[4] = pm.insert_at_vertices(cv[4], e[0]->source(), e[3]->target());
  e[5] = pm.insert_at_vertices(cv[5], e[0]->source(), e[2]->target());

  /*             
     e3  /\  e4
        /  \
        ----
       | e5 | 
    e2 |    | e0
       |    |
        ----
         e1 
  */

  std::cout << ((pm.is_valid()) ? "map valid!" : "map invalid!") << std::endl
            << std::endl;

  // Print the map:
  std::cout << std::endl << "* * * Printing map: " << std::endl << std::endl;
  CGAL::write_pm(pm, verbose_writer, std::cout);
       
  return 0;  
}
