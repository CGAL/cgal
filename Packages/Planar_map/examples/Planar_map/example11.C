// file: examples/Planar_map/example11.C

/*! \file
 * A construction of a planar map using quick insert operations. They are
 * applied through two overloaded insert methods. They are more efficient, as
 * the search for the previous halfedge in the circular list of halfedges that
 * share the same target vertex is unnecessary, saving its computation time.
 * In this example a planar map similar to the planar map of example 1 and 2
 * is constructed. The Planar_map out-streaming feature is demonstrated as
 * well. Note the format of the stream.
 */

#include "short_names.h"

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Filtered_kernel.h>
#include <CGAL/Pm_segment_traits_2.h>
#include <CGAL/Pm_default_dcel.h>
#include <CGAL/Planar_map_2.h>
#include <CGAL/IO/Pm_file_writer.h>
#include <CGAL/IO/write_pm.h>
#include <iostream>

typedef CGAL::Simple_cartesian<double>          Cartesian_kernel;
typedef CGAL::Filtered_kernel<Cartesian_kernel> Kernel;
typedef CGAL::Pm_segment_traits_2<Kernel>       Traits;
typedef Traits::Point_2                         Point_2;
typedef Traits::X_monotone_curve_2              X_monotone_curve_2;
typedef CGAL::Pm_default_dcel<Traits>           Dcel;
typedef CGAL::Planar_map_2<Dcel,Traits>         Planar_map;
typedef CGAL::Pm_file_writer<Planar_map>        Pm_writer;

int main()
{
  // Create an instance of a Planar_map:
  Planar_map pm;
  X_monotone_curve_2 cv[5];
  Planar_map::Halfedge_handle e[5];  

  Point_2 a0(100, 0), a1(20, 50), a2(100, 100), a3(180, 50);

  // Create the curves:
  cv[0] = X_monotone_curve_2(a0, a1);
  cv[1] = X_monotone_curve_2(a1, a2);
  cv[2] = X_monotone_curve_2(a2, a3);
  cv[3] = X_monotone_curve_2(a3, a0);
  cv[4] = X_monotone_curve_2(a1, a3);

  std::cout << "The curves of the map :" << std::endl;
  std::copy(&cv[0], &cv[5],
            std::ostream_iterator<X_monotone_curve_2>(std::cout, "\n"));
  std::cout << std::endl;

  // Insert the curves into the Planar_map:
  std::cout << "Inserting the curves to the map ... ";

  e[0] = pm.insert_in_face_interior(cv[0], pm.unbounded_face());
  e[1] = pm.insert_from_vertex(cv[1], e[0]);
  e[2] = pm.insert_from_vertex(cv[2], e[1]);
  e[3] = pm.insert_at_vertices(cv[3], e[2], e[0]->twin());
  e[4] = pm.insert_at_vertices(cv[4], e[1]->twin(), e[3]->twin());

  std::cout << ((pm.is_valid()) ? "map valid!" : "map invalid!") << std::endl
            << std::endl;
  
  Pm_writer verbose_writer(std::cout, pm, true);
  CGAL::write_pm(pm, verbose_writer, std::cout);

  return 0;
}
