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

#include "Map_overlay_base_test.h"

typedef CGAL::Arr_leda_segment_exact_traits  Traits;

typedef Traits::Point                                 Point;
typedef Traits::X_curve                               X_curve;
typedef Traits::Curve                                 Curve;

typedef CGAL::Map_overlay_default_dcel<Traits>        Dcel;
typedef CGAL::Planar_map_2<Dcel, Traits>              PM;

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

  if (argc < 1 || argc > 2) {
    std::cout << "usage: test data_file" << std::endl;
    exit(1);
  }

  test.start(argv[1]);
  return 0;
}

#endif // CGAL_ARR_TEST_LEDA_CONFLICT

