//#define CGAL_NOTF_DEBUG
//#define  CGAL_SWEEP_LINE_DEBUG

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
#include <CGAL/leda_real.h>
#include <CGAL/Map_overlay_default_dcel.h>
#include <CGAL/Arr_segment_circle_traits.h>
#include <CGAL/Pm_with_intersections.h>
#include <CGAL/Map_overlay_incremental.h>
#include <iostream>
#include "Map_overlay_base_test.h"

typedef leda_real                            NT;
typedef CGAL::Arr_segment_circle_traits<NT>  Traits;
typedef Traits::Segment                      Segment;
typedef Traits::Circle                       Circle;

typedef Traits::Point                                 Point;
typedef Traits::X_curve                               X_curve;
typedef Traits::Curve                                 Curve;

typedef CGAL::Map_overlay_default_dcel<Traits>        Dcel;
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
  
  Curve operator()(std::istream& is) 
  {
    Curve cv;
    
    // Get the arc type.
    char type;
    
    // Currently expects no comments in input file
    // Should be changed?
    is >> type;
    
  // A full circle (c) or a circular arc (a):
    if (type == 'c' || type == 'C' || type == 'a' || type == 'A')
      {  
        // Read the circle, using the format "x0 y0 r^2"
        NT     x0, y0, r2;
        
        is >> x0 >> y0 >> r2;
        
        Circle circle (Point (x0, y0), r2, CGAL::CLOCKWISE);
        
        if (type == 'c' || type == 'C')
          {
      // Create a full circle.
            cv = Curve(circle);  
          }
        else
          {
            // Read the end points of the circular arc.
            NT    x1, y1, x2, y2;
            
            is >> x1 >> y1 >> x2 >> y2;
            
            if ((x1 - x0)*(x1 - x0) + (y1 - y0)*(y1 - y0) != r2)
              y1 = CGAL::sqrt(r2 - (x1 - x0)*(x1 - x0)) + y0;
            
            if ((x2 - x0)*(x2 - x0) + (y2 - y0)*(y2 - y0) != r2)
              y2 = CGAL::sqrt(r2 - (x2 - x0)*(x2 - x0)) + y0;
            
            Point source (x1, y1);
            Point target (x2, y2);
            
            // Create the circular arc.
            cv = Curve (circle, source, target);
          }
      }
    else if (type == 's' || type == 'S')
      {
        // Read the end points of the segment.
        NT    x1, y1, x2, y2;
        
        is >> x1 >> y1 >> x2 >> y2;
        
        Point source (x1, y1);
        Point target (x2, y2);
        
        cv = Curve (Segment (source, target));
      }
    else
      {
        // Illegal type!
        std::cout << "Failed to read curve." << std::endl;
      }
    
    std::cout << "The read curve: " << cv << std::endl;
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

  MapOverlay_incremental              ovl_incremental;
  test.start(argv[1], &ovl_incremental);
  return 0;
}

#endif // CGAL_ARR_TEST_LEDA_CONFLICT







