// file: examples/Arrangement_2/example12.C

// Shorten long names for problematic compilers (e.g., MSVC):
#include "short_names.h"

#include <CGAL/Cartesian.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Quotient.h> 
#include <CGAL/Arr_2_default_dcel.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Arr_segment_cached_traits_2.h>
#include <CGAL/Arr_polyline_traits_2.h>
#include <iostream>
#include <CGAL/IO/Arr_iostream.h>

typedef CGAL::Quotient<CGAL::MP_Float>                  NT;
typedef CGAL::Cartesian<NT>                             Kernel;
typedef CGAL::Arr_segment_cached_traits_2<Kernel>       Seg_traits;
typedef CGAL::Arr_polyline_traits_2<Seg_traits>         Traits;

typedef Traits::Point_2                                 Point_2;
typedef Traits::Curve_2                                 Curve_2;
typedef Traits::X_monotone_curve_2                      X_monotone_curve_2;

typedef CGAL::Arr_2_default_dcel<Traits>                Dcel;
typedef CGAL::Arrangement_2<Dcel,Traits>                Arr_2;
typedef CGAL::Arr_file_writer<Arr_2>                    Arr_writer;

CGAL_BEGIN_NAMESPACE

std::ostream& operator<<(std::ostream& os, const Curve_2& cv)
{
  Curve_2::const_iterator iter;
  
  os << cv.points() << std::endl;
  for (iter = cv.begin(); iter != cv.end(); iter++)
    os << " " << *iter;

  return os;
}

std::istream& operator>>(std::istream& in, Curve_2& cv)
{
  Kernel::Point_2            p;
  std::size_t                size;
  std::list<Kernel::Point_2> pts;
  unsigned int               i;

  in >> size;

  for (i = 0; i < size; i++)
  {
    in >> p;
    pts.push_back(p);  
  }
  cv = Curve_2(pts.begin(), pts.end());
  
  return in;
}

CGAL_END_NAMESPACE

int main()
{
  Arr_2 arr;
  std::cin >> arr;
  std::cout << " * * * Printing list of all halfedges " 
            << "of the resulting Arrangement" 
            << std::endl;
  Arr_writer verbose_writer(std::cout, arr, true);
  verbose_writer.write_halfedges(arr.halfedges_begin(), arr.halfedges_end());
  return 0;
}
