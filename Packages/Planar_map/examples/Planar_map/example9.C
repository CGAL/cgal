// examples/Planar_map/example9.C
// ------------------------------

#include "short_names.h"

#include <CGAL/Cartesian.h>
#include <CGAL/Quotient.h>
#include <CGAL/Pm_default_dcel.h>
#include <CGAL/Planar_map_2.h>
#include <CGAL/Pm_segment_traits_2.h>
#include <CGAL/IO/write_pm.h>
#include <CGAL/IO/Pm_iostream.h>
#include <iostream>

// #define CGAL_POSTSCRIPT
#if defined(CGAL_USE_LEDA) && defined(CGAL_POSTSCRIPT)
#include <CGAL/IO/Pm_Postscript_file_stream.h>
#endif

typedef CGAL::Quotient<int>                     NT;
typedef CGAL::Cartesian<NT>                     Kernel;
typedef CGAL::Pm_segment_traits_2<Kernel>       Traits;
typedef CGAL::Pm_default_dcel<Traits>           Dcel;
typedef CGAL::Planar_map_2<Dcel,Traits>         Planar_map;
typedef CGAL::Pm_file_writer<Planar_map>        Pm_writer;

int main()
{ 
  Planar_map pm;
  Pm_writer verbose_writer(std::cout, pm, true);
  Pm_writer writer(std::cout, pm);

  std::cout << "* * * Demonstrating a trivial use of IO functions"
            << std::endl << std::endl;
  std::cin  >> pm;
  std::cout << pm;
  
  std::cout << std::endl;
  std::cout << "* * * Presenting the use of verbose format" << std::endl;
  std::cout << std::endl;
  CGAL::write_pm(pm, verbose_writer, std::cout);
  
  std::cout << std::endl;
  std::cout << "* * * Demonstrating the use of the writer class interface."
            << std::endl;
  std::cout << "* * * Printing all halfedges in non verbose format"
            << std::endl << std::endl;
  writer.write_halfedges(pm.halfedges_begin(), pm.halfedges_end());
  std::cout << std::endl;
  std::cout << "* * * Printing all halfedges in a verbose format" << std::endl
            << std::endl;
  verbose_writer.write_halfedges(pm.halfedges_begin(), pm.halfedges_end());
   
#if defined(CGAL_USE_LEDA) && defined(CGAL_POSTSCRIPT)
  // Print to Postscript file:
  CGAL::Postscript_file_stream  LPF(500, 500 ,"pm.ps");
  LPF.init(-3,3,-3);
  LPF.set_line_width(1);
  LPF << pm;
#endif

  return 0;
}
