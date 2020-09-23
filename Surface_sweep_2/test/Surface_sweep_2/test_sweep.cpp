// examples/Pm_with_intersections/example4
// ---------------------------------------
#include <CGAL/config.h>
#include <algorithm>
#include <list>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>

#define CGAL_SEGMENT_TRAITS          1
#define CGAL_SEGMENT_LEDA_TRAITS     2
#define CGAL_POLYLINE_TRAITS         11
#define CGAL_CONIC_TRAITS            21

// Picking a default Traits class (this, with the
// PL flag enables the running of the test independently of cgal_make.)
#ifndef CGAL_ARR_TEST_TRAITS
#define CGAL_ARR_TEST_TRAITS CGAL_SEGMENT_TRAITS
#endif

// Making sure test doesn't fail if LEDA is not installed
#if ! defined(CGAL_USE_LEDA) && \
      (CGAL_ARR_TEST_TRAITS == CGAL_SEGMENT_LEDA_TRAITS || \
       CGAL_ARR_TEST_TRAITS == CGAL_CONIC_TRAITS )

int main()
{
  std::cout << "A try to run test with LEDA traits but LEDA is not installed.";
  std::cout << std::endl;
  std::cout << "Test is not performed.";
  std::cout << std::endl;

  return 0;
}
#elif ! defined(CGAL_USE_GMP) && \
        (CGAL_ARR_TEST_TRAITS == CGAL_SEGMENT_TRAITS)

int main()
{
  std::cout << "A try to run test with GMP number type but GMP is not installed.";
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
#include<CGAL/Gmpq.h>
#elif CGAL_ARR_TEST_TRAITS == CGAL_SEGMENT_LEDA_TRAITS
#include <CGAL/leda_rational.h>
#include <CGAL/Pm_segment_traits_leda_kernel_2.h>
#include <CGAL/Arr_leda_segment_traits_2.h>
#elif CGAL_ARR_TEST_TRAITS == CGAL_POLYLINE_TRAITS
#include <CGAL/Cartesian.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Quotient.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arr_polyline_traits_2.h>
#elif CGAL_ARR_TEST_TRAITS == CGAL_CONIC_TRAITS
#include <CGAL/Cartesian.h>
#include <CGAL/leda_real.h>
#include <CGAL/Arr_conic_traits_2.h>
#else
#error No traits defined for test
#endif

#include <list>
#include <string>
#include <CGAL/Surface_sweep_2_algorithms.h>
#include "CompareCurveList.h"

#if CGAL_ARR_TEST_TRAITS==CGAL_SEGMENT_TRAITS
  typedef CGAL::Gmpq                                    NT;
  typedef CGAL::Cartesian<NT>                           Kernel;
  typedef CGAL::Arr_segment_traits_2<Kernel>            Traits;

#elif CGAL_ARR_TEST_TRAITS == CGAL_SEGMENT_LEDA_TRAITS
  typedef leda_rational                                 NT;
  typedef CGAL::Pm_segment_traits_leda_kernel_2         Kernel;
  typedef CGAL::Arr_leda_segment_traits_2<Kernel>       Traits;

#elif CGAL_ARR_TEST_TRAITS == CGAL_POLYLINE_TRAITS
  typedef CGAL::Quotient<MP_Float>                      NT;
  typedef CGAL::Cartesian<NT>                           Kernel;
  typedef CGAL::Arr_segment_cached_traits_2<Kernel>     Seg_traits;
  typedef CGAL::Arr_polyline_traits_2<Seg_traits>       Traits;

#elif CGAL_ARR_TEST_TRAITS == CGAL_CONIC_TRAITS
typedef leda_real                                       NT;
typedef CGAL::Cartesian<NT>                             Kernel;
typedef CGAL::Arr_conic_traits_2<Kernel>                Traits;
#endif

typedef Traits::Point_2                     Point_2;
typedef Traits::X_monotone_curve_2                   X_monotone_curve_2;

typedef std::list<Point_2> PointList;
typedef PointList::iterator PointListIter;

typedef std::list<X_monotone_curve_2> CurveList;
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

  Traits tr;

  // get subcurves w/o overlapping
  CurveList curves_no_overlap_list_out;
  CGAL::compute_subcurves(curves.begin(),
                          curves.end(),
                          std::back_inserter(curves_no_overlap_list_out));


  // get subcurves w/ overlapping
  CurveList curves_with_overlap_list_out;
  CGAL::compute_subcurves(curves.begin(),
                          curves.end(),
                          std::back_inserter(curves_with_overlap_list_out),
                          true);

  /*std::copy(curves_no_overlap_list_out.begin(),
            curves_no_overlap_list_out.end(),
            std::ostream_iterator<X_monotone_curve_2>(std::cout, "\n"));
  std::cout<<"\n\n*******************\n\n";

  std::copy(curves_with_overlap_list_out.begin(),
            curves_with_overlap_list_out.end(),
            std::ostream_iterator<X_monotone_curve_2>(std::cout, "\n"));
  return 0;*/

  //std::cout << mylist1.size() << " curves\n";


  // get intersection points (with endpoints)
  PointList points_with_ends_list_out;
  CGAL::compute_intersection_points(curves.begin(),
                                    curves.end(),
                                    std::back_inserter(points_with_ends_list_out),
                                    true);


  // get intersection points w/o end points
  PointList points_without_ends_list_out;
  CGAL::compute_intersection_points(curves.begin(),
                                    curves.end(),
                                    std::back_inserter(points_without_ends_list_out),
                                    false);
  std::cout << points_without_ends_list_out.size()
            << " points_without_ends_list_out(size)\n";

  // check the do_curves_intersecting method
  bool do_intersect_out =
    CGAL::do_curves_intersect(curves.begin(), curves.end());


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

  if ( !compare_lists(curves_no_overlap_list_out,
                      curves_no_overlap_list, tr) )
    return -1;

  if ( !compare_lists(curves_with_overlap_list_out,
                      curves_with_overlap_list, tr) )
    return -1;

  if ( !compare_lists(points_with_ends_list_out,
                      points_with_ends_list, tr))
    return -1;

  if ( !compare_lists(points_without_ends_list_out,
                      points_without_ends_list, tr))
    return -1;

  bool do_intersect = false;
  if((points_without_ends_list.size() != 0) ||
     (curves_no_overlap_list_out.size() !=
      curves_with_overlap_list_out.size()))
    do_intersect = true;

  if (do_intersect_out != do_intersect)
    return -1;

  std::cout<<"OK\n";
  return 0;
}

void ReadCurveList(std::ifstream& inp, CurveList& clist)
{
  int count;
  inp >> count;
  //std::cout << "ReadCurveList " << count << "\n";

  for (int i = 0; i < count; i++) {
    NT x0, y0, x1, y1;
    int ix0, iy0, ix1, iy1;
    inp >> ix0 >> iy0 >> ix1 >> iy1;
    x0 = ix0; y0 = iy0; x1 = ix1; y1 = iy1;

    Point_2 p1(x0, y0);
    Point_2 p2(x1, y1);
    X_monotone_curve_2 curve(p1, p2);
    clist.push_back(curve);
    //std::cout << curve << "\n";
  }
}

void ReadCurveListRational(std::ifstream& inp, CurveList& clist)
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

    X_monotone_curve_2 curve(p1, p2);
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

#endif
