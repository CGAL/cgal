// see examples/Planar_map/example1.C
// for original version ...
// ------------------------------

// needed in Pm_segment_traits_2 - otherwise
// we use x() / y()
#define NBUG
#define LEDA_NO_MIN_MAX_TEMPL
#define CGAL_PROVIDE_LEDA_RAT_KERNEL_TRAITS_3

#include <CEP/Leda_rat_kernel/leda_rat_kernel_traits.h>
#include <CGAL/Kernel_checker.h>
#include <CGAL/Pm_segment_traits_2.h>
#include <CGAL/Pm_default_dcel.h>
#include <CGAL/Planar_map_2.h>
#include <CGAL/IO/Pm_file_writer.h>
#include <CGAL/IO/write_pm.h>
#include <iostream>

//typedef CGAL::leda_rat_kernel_traits      K1;
//typedef CGAL::Homogeneous<leda_integer>   K2;
//typedef CGAL::Kernel_checker<K1, K2, CGAL::leda_to_cgal_2 > Kernel;

typedef CGAL::leda_rat_kernel_traits            Kernel;
typedef CGAL::Pm_segment_traits_2<Kernel>       Traits;
typedef Traits::Point_2                         Point_2;
typedef Traits::X_curve_2                       X_curve_2;
typedef CGAL::Pm_default_dcel<Traits>           Dcel;
typedef CGAL::Planar_map_2<Dcel,Traits>         Planar_map;
typedef CGAL::Pm_file_writer<Planar_map>        Pm_writer;

int main()
{
  // Create an instance of a Planar_map:
  Planar_map pm;
  Pm_writer verbose_writer(std::cout, pm, true);
  X_curve_2 cv[5];
  int i;

  CGAL::set_ascii_mode(std::cout);

  Point_2 a1(100, 0), a2(20, 50), a3(180, 50), a4(100, 100);

  // Create the curves:
  cv[0] = X_curve_2(a1, a2);
  cv[1] = X_curve_2(a1, a3);
  cv[2] = X_curve_2(a2, a3);
  cv[3] = X_curve_2(a2, a4);
  cv[4] = X_curve_2(a3, a4);
  
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
  X_curve_2 c1(a2, p);
  X_curve_2 c2(p, a3);
  Planar_map::Halfedge_handle se = pm.split_edge(e[2], c1, c2); 
  pm.insert_at_vertices(X_curve_2(p, a1), se->target(), e[0]->source());

  // Print map after splitting and adding:
  std::cout << std::endl << "* * * Map after:" << std::endl << std::endl;
  CGAL::write_pm(pm, verbose_writer, std::cout);

  return 0;  
}
