//#define CGAL_NOTF_DEBUG
#include <CGAL/config.h> // needed for the LONGNAME flag

#ifdef CGAL_CFG_NO_LONGNAME_PROBLEM
// Define shorter names to please linker (g++/egcs)
#define Arrangement_2 Ar
#define Cartesian Cr
#define Arr_polyline_traits APT
#define Quotient Qt
#define Planar_ma2p_traits_wrap PMTW
#define _List_iterator Lit
#define bidirectional_iterator_tag Bitt
#define Planar_map_2 Pm2
#define Arr_2_default_dcel A2dd
#define Point_2 Pt2
#define allocator Altr
#define Td_X_trapezoid TdXt
#define Td_traits Tdt
#endif

#include <CGAL/Cartesian.h>
//#include <CGAL/Arr_2_bases.h>
//#include <CGAL/Arr_2_default_dcel.h>
#include <fstream>

#ifndef CGAL_PLANAR_MAP_2
#include <CGAL/Planar_map_2.h>
#endif

#ifndef CGAL_ARR_2_OVERLAY_DCEL_H
#include <CGAL/Map_overlay_default_dcel.h>
#endif

#ifndef CGAL_MAP_OVERLAY_POST_PROC_NOTIFIER_H
#include <CGAL/Map_overlay_post_proc_notifier.h>
#endif

#ifndef CGAL_MAP_OVERLAY_H
#include <CGAL/Map_overlay.h>
#endif

// Picking a default Traits class (this, with the 
// PL flag enables the running of the test independently of cgal_make.)
#ifndef CGAL_ARR_TEST_TRAITS
//#define CGAL_ARR_TEST_TRAITS CGAL_SEGMENT_TRAITS
#define CGAL_ARR_TEST_TRAITS CGAL_SEGMENT_LEDA_TRAITS
//#define CGAL_ARR_TEST_TRAITS CGAL_POLYLINE_TRAITS
//#define CGAL_ARR_TEST_TRAITS CGAL_POLYLINE_LEDA_TRAITS
//#define CGAL_ARR_TEST_TRAITS CGAL_SEGMENT_CIRCLE_TRAITS
#endif

// Making sure test doesn't fail if LEDA is not installed
#if ! defined(CGAL_USE_LEDA) && \
      (CGAL_ARR_TEST_TRAITS == CGAL_POLYLINE_LEDA_TRAITS || \
       CGAL_ARR_TEST_TRAITS == CGAL_SEGMENT_LEDA_TRAITS || \
       CGAL_ARR_TEST_TRAITS == CGAL_SEGMENT_CIRCLE_TRAITS)

int main(int argc, char* argv[])
{
  std::cout << "A try to run test with LEDA traits but LEDA is not installed.";
  std::cout << std::endl;
  std::cout << "Test is not performed.";
  std::cout << std::endl;

  return 0;
}
#else

// Choose traits

#include <CGAL/leda_rational.h>
#include <CGAL/Arr_leda_segment_exact_traits.h>

#include <CGAL/Pm_walk_along_line_point_location.h>

 
//#include <CGAL/Arrangement_2.h>

#include "Map_overlay_base_test.h"

// Quotient is included anyway, because it is used to read
// data files. Quotient can read both integers and fractions.
// leda rational will only read fractions.
#include <CGAL/Quotient.h> 

#include <list>
#include <string>


typedef leda_rational                        NT;
typedef CGAL::Arr_leda_segment_exact_traits  Traits;

typedef Traits::Point                                 Point;
typedef Traits::X_curve                               X_curve;
typedef Traits::Curve                                 Curve;

typedef CGAL::Map_overlay_default_dcel<Traits>                Dcel;

typedef CGAL::Planar_map_2<Dcel, Traits>                      PM;

class Read_segment {
public:
  Curve operator()(std::istream& is) {
    int   x1,y1,x2,y2;
    
    is >> x1 >> y1 >> x2 >> y2;
    
    Point p1(x1,y1), p2(x2,y2);
    Curve cv(p1,p2);
    
    return cv;
  }
};


int main(int argc, char* argv[])
{
  
  Map_overlay_base_test<PM, Read_segment> test;

  if (argc < 2 || argc > 3) {
    std::cout << "usage: test data_file" << std::endl;
    exit(1);
  }

  test.start(argv[1], argv[2]);
  return 0;
}

#endif // CGAL_ARR_TEST_LEDA_CONFLICT
