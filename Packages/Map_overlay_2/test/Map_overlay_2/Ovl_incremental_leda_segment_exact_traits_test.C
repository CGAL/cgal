//#define  CGAL_SWEEP_LINE_DEBUG
//#define  OVL_DEBUG_TEST

#include <CGAL/config.h> // needed for the LONGNAME flag
#include <iostream>

// Making sure test doesn't fail if LEDA is not installed
#if ! defined(CGAL_USE_LEDA)
int main(int argc, char* argv[])
{
  std::cout << "A try to run test with LEDA traits but LEDA is not installed.";
  std::cout << std::endl;
  std::cout << "Test is not performed.";
  std::cout << std::endl;

  return 0;
}
#else

//#include <CGAL/Cartesian.h>
#include <CGAL/Arr_leda_segment_exact_traits.h>

#include <CGAL/Pm_with_intersections.h>

#include "Map_overlay_base_test.h"

#ifndef CGAL_MAP_OVERLAY_INCREMENTAL_H
#include <CGAL/Map_overlay_incremental.h>
#endif

typedef CGAL::Arr_leda_segment_exact_traits  Traits;

typedef Traits::Point                                 Point;
typedef Traits::X_curve                               X_curve;
typedef Traits::Curve                                 Curve;

typedef CGAL::Map_overlay_default_dcel<Traits>        Dcel;
typedef CGAL::Planar_map_2<Dcel, Traits>              Planar_map;
typedef CGAL::Planar_map_with_intersections_2<Planar_map>   Pmwx;
typedef CGAL::Map_overlay_default_notifier<Planar_map>      
                                                 MapOverlay_change_notification;
typedef CGAL::Map_overlay_incremental<Pmwx,MapOverlay_change_notification>   
                                                               MapOverlay_incremental;

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
  
  Map_overlay_base_test<Pmwx, Read_segment> test;

  if (argc < 1 || argc > 2) {
    std::cout << "usage: test data_file" << std::endl;
    exit(1);
  }

  MapOverlay_incremental    ovl_incremental;
  test.start(argv[1], &ovl_incremental);
  return 0;
}

#endif // CGAL_ARR_TEST_LEDA_CONFLICT





