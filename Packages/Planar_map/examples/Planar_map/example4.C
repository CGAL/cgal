/*! \file
 * A construction of a planar map using various insert and split operations.
 * In this example a planar map similar to the planar map of example 1 and 2
 * is constructed. The planar map is further modified by spliting one of its
 * edges and inserting an additional edge at specified vertices. The
 * Planar_map out-streaming feature is demonstrated as well. Note the format
 * of the stream.
 */

#include "short_names.h"

#include <CGAL/Homogeneous.h>
#include <CGAL/Pm_segment_traits_2.h>
#include <CGAL/Pm_default_dcel.h>
#include <CGAL/Planar_map_2.h>
#include <CGAL/IO/Pm_file_writer.h>
#include <CGAL/IO/write_pm.h>
#include <iostream>

typedef CGAL::Homogeneous<long>                 Kernel;
typedef CGAL::Pm_segment_traits_2<Kernel>       Traits;
typedef Traits::Point_2                         Point_2;
typedef Traits::X_monotone_curve_2                       X_monotone_curve_2;
typedef CGAL::Pm_default_dcel<Traits>           Dcel;
typedef CGAL::Planar_map_2<Dcel,Traits>         Planar_map;
typedef CGAL::Pm_file_writer<Planar_map>        Pm_writer;

int main()
{
  // Create an instance of a Planar_map:
  Planar_map pm;
  Pm_writer verbose_writer(std::cout, pm, true);
  X_monotone_curve_2 cv[5];
  int i;

  CGAL::set_ascii_mode(std::cout);

  Point_2 a1(100, 0), a2(20, 50), a3(180, 50), a4(100, 100);

  // Create the curves:
  cv[0] = X_monotone_curve_2(a1, a2);
  cv[1] = X_monotone_curve_2(a1, a3);
  cv[2] = X_monotone_curve_2(a2, a3);
  cv[3] = X_monotone_curve_2(a2, a4);
  cv[4] = X_monotone_curve_2(a3, a4);
  
  // Insert the curves into the Planar_map:
  std::cout << "Inserting the curves to the map ... ";
  Planar_map::Halfedge_handle e[5];  
  for (i = 0; i < 5; i++)
    e[i] = pm.insert(cv[i]);
  std::cout << ((pm.is_valid()) ? "map valid!" : "map invalid!") << std::endl
            << std::endl;
  
  // Print map before splitting and adding:
  std::cout << "* * * Map before:" << std::endl << std::endl;
  CGAL::write_pm(pm, verbose_writer, std::cout);
  
  // Split e[2] in the middle, and add a curve between the new vertex and
  // the source of e[0]:
  Point_2 p(100, 50);
  X_monotone_curve_2 c1(a2, p);
  X_monotone_curve_2 c2(p, a3);
  Planar_map::Halfedge_handle se = pm.split_edge(e[2], c1, c2); 
  pm.insert_at_vertices(X_monotone_curve_2(p, a1), se->target(),
                        e[0]->source());

  // Print map after splitting and adding:
  std::cout << std::endl << "* * * Map after:" << std::endl << std::endl;
  CGAL::write_pm(pm, verbose_writer, std::cout);

  return 0;  
}
