// file: examples/Arrangement_2/example11.C

#include "short_names.h"

#include <CGAL/Cartesian.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Quotient.h>
#include <CGAL/Arr_2_default_dcel.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/IO/Arr_iostream.h>
#include <iostream>

#ifdef CGAL_USE_LEDA
// #include <CGAL/IO/Arr_Postscript_file_stream.h>
#endif

typedef CGAL::Quotient<CGAL::MP_Float>                  NT;
typedef CGAL::Cartesian<NT>                             Kernel;
typedef CGAL::Arr_segment_traits_2<Kernel>              Traits;
typedef Traits::Curve                                   Curve;
typedef Traits::X_monotone_curve_2                      X_monotone_curve_2;
typedef CGAL::Arr_2_default_dcel<Traits>                Dcel;
typedef CGAL::Arrangement_2<Dcel,Traits>                Arr;

int main()
{
  Arr arr;

  std::cout << "* * * Demonstrating a trivial use of IO functions";
  std::cout << std::endl << std::endl;
  std::cin >> arr;
  std::cout << arr;
  
  std::cout << std::endl;

  std::cout << "* * * Presenting the use of verbose format";
  std::cout << std::endl << std::endl;;
  CGAL::Arr_file_writer<Arr> verbose_writer(std::cout, arr, true);
  CGAL::write_arr(arr, verbose_writer, std::cout);

  // printing to Postscript file.
#ifdef CGAL_USE_LEDA
  //  CGAL::Postscript_file_stream  LPF(500, 500 ,"arr.ps");
  //  LPF.init(-3,3,-3);
  //  LPF.set_line_width( 1);
  //  LPF << arr;
#endif

  return 0;
}
