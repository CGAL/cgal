// examples/Pm_with_intersections/example4
// ---------------------------------------
#include <algorithm>
#include <list>
#include <iostream>
#include <sstream>
#include <string>

#include "short_names.h"

#include <CGAL/basic.h>

#define CGAL_SEGMENT_TRAITS          1
#define CGAL_SEGMENT_LEDA_TRAITS     2
#define CGAL_POLYLINE_TRAITS         11
#define CGAL_POLYLINE_LEDA_TRAITS    12
#define CGAL_SEGMENT_CIRCLE_TRAITS   21

// Picking a default Traits class (this, with the 
// PL flag enables the running of the test independently of cgal_make.)
#ifndef CGAL_ARR_TEST_TRAITS
#define CGAL_ARR_TEST_TRAITS CGAL_SEGMENT_LEDA_TRAITS
#endif

// Making sure test doesn't fail if LEDA is not installed
#if ! defined(CGAL_USE_LEDA) && \
      (CGAL_ARR_TEST_TRAITS == CGAL_POLYLINE_LEDA_TRAITS || \
       CGAL_ARR_TEST_TRAITS == CGAL_SEGMENT_LEDA_TRAITS || \
       CGAL_ARR_TEST_TRAITS == CGAL_SEGMENT_CIRCLE_TRAITS )

int main()
{
  std::cout << "A try to run test with LEDA traits but LEDA is not installed.";
  std::cout << std::endl;
  std::cout << "Test is not performed.";
  std::cout << std::endl;

  return 0;
}
#else

// Choose traits

#if CGAL_ARR_TEST_TRAITS==CGAL_SEGMENT_TRAITS 
#include <CGAL/Cartesian.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Quotient.h>
#include <CGAL/Arr_segment_traits_2.h>
#elif CGAL_ARR_TEST_TRAITS == CGAL_SEGMENT_LEDA_TRAITS
#include <CGAL/leda_rational.h>
#include <CGAL/Pm_segment_traits_leda_kernel_2.h>
#include <CGAL/Arr_leda_segment_traits_2.h>
#elif CGAL_ARR_TEST_TRAITS == CGAL_POLYLINE_TRAITS
#include <CGAL/Cartesian.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Quotient.h>
#include <CGAL/Arr_polyline_traits.h>
#elif CGAL_ARR_TEST_TRAITS == CGAL_POLYLINE_LEDA_TRAITS
#include <CGAL/leda_rational.h>
#include <CGAL/Pm_segment_traits_leda_kernel_2.h>
#include <CGAL/Arr_leda_polyline_traits.h>
#elif CGAL_ARR_TEST_TRAITS == CGAL_SEGMENT_CIRCLE_TRAITS
#include <CGAL/leda_real.h>
#include <CGAL/Arr_segment_circle_traits.h>
#else
#error No traits defined for test
#endif

#include <list>
#include <string>
#include <CGAL/Sweep_line_tight_2.h>
#include <CGAL/Sweep_line_2/Sweep_line_event.h>
#include <CGAL/Sweep_line_2/Sweep_line_subcurve.h>
#include "CompareCurveList.h"

#if CGAL_ARR_TEST_TRAITS==CGAL_SEGMENT_TRAITS 
  typedef CGAL::Quotient<int>                           NT;
  typedef CGAL::Cartesian<NT>                           Kernel;
  typedef CGAL::Arr_segment_traits_2<Kernel>            Traits;

#elif CGAL_ARR_TEST_TRAITS == CGAL_SEGMENT_LEDA_TRAITS
  typedef leda_rational                                 NT;
  typedef CGAL::Pm_segment_traits_leda_kernel_2         Kernel;
  typedef CGAL::Arr_leda_segment_traits_2<Kernel>       Traits;

#elif CGAL_ARR_TEST_TRAITS == CGAL_POLYLINE_TRAITS
  typedef CGAL::Quotient<int>                           NT;
  typedef CGAL::Cartesian<NT>                           Kernel;
  typedef CGAL::Arr_polyline_traits<Kernel>             Traits;

#elif CGAL_ARR_TEST_TRAITS == CGAL_POLYLINE_LEDA_TRAITS
  typedef leda_rational                                 NT;
  typedef CGAL::Pm_segment_traits_leda_kernel_2         Kernel;
  typedef CGAL::Arr_leda_polyline_traits<Kernel>        Traits;

#elif CGAL_ARR_TEST_TRAITS == CGAL_SEGMENT_CIRCLE_TRAITS
  typedef leda_real                                     NT;
  typedef CGAL::Arr_segment_circle_traits<NT>           Traits;
  typedef Traits::Segment                               Segment;
  typedef Traits::Circle                                Circle;
#endif

typedef Traits::Point_2                     Point_2;
typedef Traits::X_curve_2                   X_curve_2;

typedef std::list<Point_2> PointList;
typedef PointList::iterator PointListIter;

typedef CGAL::Sweep_line_subcurve<Traits> SubCurve;
typedef CGAL::Sweep_line_event<Traits, SubCurve> Event;

typedef std::list<X_curve_2> CurveList;
typedef CurveList::iterator CurveListIter;

void ReadCurveList(std::ifstream &inp, CurveList &clist);
void ReadCurveListRational(std::ifstream &inp, CurveList &clist);
void ReadPointList(std::ifstream &inp, PointList &plist);
bool IsCurveListIdentical(CurveList &list1, CurveList &list2);
bool IsPointListIdentical(PointList &list1, PointList &list2);

int main(int argc, char * argv[])
{

  if ( argc != 2 )
  {
    std::cout << "Specify a file name " << std::endl;
    return -1;
  }

  std::ifstream inp(argv[1]);
  if (!inp.is_open()) {
    std::cerr << "Cannot open file " << argv[1] << "!" << std::endl;
    return -1;
  }
  
  CurveList curves;
  ReadCurveList(inp, curves);
  
  // get subcurves w/o overlapping
  CurveList curves_no_overlap_list_out;
  CGAL::Sweep_line_tight_2<CurveListIter, Traits, Event, SubCurve> sl;
  sl.get_subcurves(curves.begin(), curves.end(), 
		   std::back_inserter(curves_no_overlap_list_out));
  //std::cout << mylist.size() << " curves\n";


  // get subcurves w/ overlapping
  CurveList curves_with_overlap_list_out;
  CGAL::Sweep_line_tight_2<CurveListIter, Traits, Event, SubCurve> sl1;
  sl1.get_subcurves(curves.begin(), curves.end(), 
		   std::back_inserter(curves_with_overlap_list_out), true);
  //std::cout << mylist1.size() << " curves\n";


  // get intersection points (with endpoints)
  PointList points_with_ends_list_out;
  CGAL::Sweep_line_tight_2<CurveListIter, Traits, Event, SubCurve> sl2;
  sl2.get_intersection_points(curves.begin(), curves.end(), 
			     std::back_inserter(points_with_ends_list_out));
  //std::cout << mypointlist.size() << " points \n";


  // get intersection points w/o end points
  PointList points_without_ends_list_out;
  CGAL::Sweep_line_tight_2<CurveListIter, Traits, Event, SubCurve> sl3;
  sl3.get_intersection_points(curves.begin(), curves.end(), 
			     std::back_inserter(points_without_ends_list_out), false);
  //std::cout << mypointlist2.size() << " points (internal)\n";


  // read curves and points from file
  CurveList curves_no_overlap_list;
  ReadCurveListRational(inp, curves_no_overlap_list);

  PointList points_with_ends_list;
  ReadPointList(inp, points_with_ends_list);

  PointList points_without_ends_list;
  ReadPointList(inp, points_without_ends_list);

  int num_faces;
  inp >> num_faces;

  CurveList curves_with_overlap_list;
  ReadCurveListRational(inp, curves_with_overlap_list);

  if ( !IsCurveListIdentical(curves_no_overlap_list_out, curves_no_overlap_list) )
    return -1;

  if ( !IsCurveListIdentical(curves_with_overlap_list_out, curves_with_overlap_list) )
    return -1;

  if ( !IsPointListIdentical(points_with_ends_list_out, points_with_ends_list))
    return -1;

  if ( !IsPointListIdentical(points_without_ends_list_out, points_without_ends_list))
    return -1;

  return 0;
}

void ReadCurveList(std::ifstream &inp, CurveList &clist)
{
  int count;
  inp >> count;
  std::cout << "ReadCurveList " << count << "\n";

  for (int i = 0; i < count; i++) {
    NT x0, y0, x1, y1;
    int ix0, iy0, ix1, iy1;
    inp >> ix0 >> iy0 >> ix1 >> iy1;
    x0 = ix0; y0 = iy0; x1 = ix1; y1 = iy1;

    Point_2 p1(x0, y0);
    Point_2 p2(x1, y1);
    X_curve_2 curve(p1, p2);
    clist.push_back(curve);
    std::cout << curve << "\n";
  }
}

void ReadCurveListRational(std::ifstream &inp, CurveList &clist)
{
  int count;
  inp >> count;
  std::cout << "ReadCurveListRational " << count << "\n";
  char ch;

  for (int i = 0; i < count; i++) {
    int a, b;
    inp >> a >> ch >> b;
    NT x0(a,b);
    inp >> a >> ch >> b;
    NT y0(a,b);
    Point_2 p1(x0, y0);

    inp >> a >> ch >> b;
    NT x1(a,b);
    inp >> a >> ch >> b;
    NT y1(a,b);
    Point_2 p2(x1, y1);

    X_curve_2 curve(p1, p2);
    clist.push_back(curve);
    std::cout << curve << "\n";
  }
}
void ReadPointList(std::ifstream &inp, PointList &plist)
{
  int count;
  inp >> count;
  char ch;

  std::cout << "ReadPointList " << count << "\n";
  for (int i = 0; i < count; i++) {
    int a, b;
    inp >> a >> ch >> b;
    NT x0(a,b);
    inp >> a >> ch >> b;
    NT y0(a,b);
    Point_2 p(x0, y0);
    plist.push_back(p);
  }
}

bool IsCurveListIdentical(CurveList &list1, CurveList &list2)
{
  CompareCurveList<CurveList, Traits> comparer;
  return comparer.IsIdentical(list1.begin(), list1.end(), list2.begin(), list2.end());
}

bool IsPointListIdentical(PointList &list1, PointList &list2)
{
  ComparePointList<PointList, Traits> comparer;
  return comparer.IsIdentical(list1.begin(), list1.end(), list2.begin(), list2.end());
}


#endif
