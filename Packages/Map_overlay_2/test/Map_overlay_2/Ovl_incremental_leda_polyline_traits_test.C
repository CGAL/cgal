//#define  CGAL_SWEEP_LINE_DEBUG
//#define  OVL_DEBUG_TEST

#include <CGAL/config.h> // needed for the LONGNAME flag

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

#include <CGAL/basic.h>
#include <CGAL/leda_rational.h> 
#include <CGAL/Map_overlay_default_dcel.h>
#include <CGAL/Arr_leda_polyline_traits.h>
#include <CGAL/Pm_with_intersections.h>
#include <CGAL/Map_overlay_incremental.h>
#include <iostream>
#include "Map_overlay_base_test.h"

typedef leda_rational                                 NT;
typedef CGAL::Arr_leda_polyline_traits<NT>            Traits;

typedef Traits::Point                                 Point;
typedef Traits::X_curve                               X_curve;
typedef Traits::Curve                                 Curve;

typedef CGAL::Map_overlay_default_dcel<Traits>         Dcel;
typedef CGAL::Planar_map_2<Dcel, Traits>               Planar_map;
typedef CGAL::Planar_map_with_intersections_2<Planar_map>   Pmwx;
typedef CGAL::Map_overlay_default_notifier<Planar_map>      
                                         MapOverlay_change_notification;
typedef CGAL::Map_overlay_incremental<Pmwx, MapOverlay_change_notification>   
                                                       MapOverlay_incremental;

// This is provided dute to the fact that the extractor (is >> ) 
// is not defined per each kind of curve. 
class Read_segment {
public:
  Curve operator()(std::istream& is) {
    Curve cv;
    
    std::size_t  size;
    
    is >> size;

    for (unsigned int i = 0; i < size; ++i){
      int x,y;
      
      is>>x>>y;
      
      cv.push_back(Point(x,y));  
    }
    return cv;
  }
};

#ifdef OVL_DEBUG_TEST
std::ostream&  operator<<(std::ostream& os, const Curve & cv)
{
  typedef Curve::const_iterator       Points_iterator;
  
  os << cv.size() << std::endl;
  
  for (Points_iterator points_iter = cv.begin(); 
       points_iter != cv.end(); points_iter++)
    os<<" "<<*points_iter;
  
  return os;
}
#endif


int main(int argc, char* argv[])
{
  
  Map_overlay_base_test<Pmwx, Read_segment> test;

  if (argc < 1 || argc > 2) {
    std::cout << "usage: test data_file" << std::endl;
    exit(1);
  }

  MapOverlay_incremental        ovl_incremental;
  test.start(argv[1], &ovl_incremental);
  return 0;
}

#endif // CGAL_ARR_TEST_LEDA_CONFLICT
