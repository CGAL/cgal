//examples/Arrangement_2/example12

#if defined(CGAL_CFG_NO_LONGNAME_PROBLEM) || defined(_MSC_VER)
// Define shorter names to please linker (g++/egcs)
#define Arrangement_2                           _Ar
#define Cartesian                               _Cr
#define Quotient                                _Qt
#define Planar_map_2                            _PM
#define allocator                               _All
#define Pm_default_dcel                         _PDD
#define Arr_polyline_traits                     _APT
#define Arr_2_default_dcel                      _ADD
#endif

#include <CGAL/Cartesian.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Quotient.h> 
#include <CGAL/Arr_2_bases.h>
#include <CGAL/Arr_2_default_dcel.h>
#include <CGAL/Arrangement_2.h>
#include <iostream>
#include <CGAL/Arr_polyline_traits.h>
#include <CGAL/IO/Arr_iostream.h>

typedef CGAL::Quotient<CGAL::MP_Float>                  NT;
typedef CGAL::Cartesian<NT>                             Kernel;
typedef CGAL::Arr_polyline_traits<Kernel>               Traits;

typedef Traits::Point                                   Point;
typedef Traits::Curve_2                                 Curve_2;

typedef CGAL::Arr_2_default_dcel<Traits>                Dcel;
typedef CGAL::Arr_base_node<Curve_2>                    Base_node;
typedef CGAL::Arrangement_2<Dcel,Traits,Base_node>      Arr_2;

typedef CGAL::Arr_file_writer<Arr_2>                    Arr_writer;

CGAL_BEGIN_NAMESPACE

std::ostream & operator<<(std::ostream & os, const Curve_2 & cv)
{
  typedef Curve_2::const_iterator Points_iterator;
  
  os << cv.size() << std::endl;
  for (Points_iterator points_iter = cv.begin(); 
       points_iter != cv.end(); points_iter++)
    os << " " << *points_iter;

  return os;
}

std::istream & operator>>(std::istream & in, Curve_2 & cv)
{
  std::size_t  size;
  unsigned int i;

  in >> size;

  for (i = 0; i < size; i++){
    Point p;
    in >> p;
    cv.push_back(p);  
  }
  
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
