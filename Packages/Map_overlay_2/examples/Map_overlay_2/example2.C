#include <CGAL/basic.h> //CGAL definitions that need to be before anything else
#include <CGAL/Cartesian.h>
#include <CGAL/Quotient.h>
#include <CGAL/Pm_walk_along_line_point_location.h>
#include <CGAL/Map_overlay_default_dcel.h>
#include <CGAL/Map_overlay.h>
#include <CGAL/Map_overlay_incremental.h>
#include <CGAL/Map_overlay_default_notifier.h>
#include <CGAL/Pm_with_intersections.h>
#include <CGAL/Arr_segment_exact_traits.h>
#include <CGAL/IO/Pm_Window_stream.h>

#include <iostream.h>
#include <vector>
#include <list>
#include <string>

#include <CGAL/IO/cgal_window.h>  //used for visualization -

typedef CGAL::Quotient<int>                    NT;
typedef CGAL::Cartesian<NT>                    R;
typedef CGAL::Arr_segment_exact_traits<R>      Traits;
typedef Traits::Point                          Point;
typedef Traits::X_curve                        X_curve;
typedef Traits::Curve                          Curve;

typedef CGAL::Map_overlay_default_dcel<Traits> Dcel;
typedef CGAL::Planar_map_2<Dcel,Traits>        Planar_map;
typedef CGAL::Planar_map_with_intersections_2<Planar_map>    Pmwx;
typedef CGAL::Map_overlay_incremental<Pmwx>     MapOverlay_incremental;
typedef CGAL::Map_overlay_2<Pmwx>               MapOverlay;

typedef CGAL::Pm_walk_along_line_point_location<Planar_map>  PmWalkPL;


using namespace std;


int  main()
{
  Pmwx  pmwx1, pmwx2;
  int   num_curves1, num_curves2;
  
  NT        x1, y1, x2, y2;

  std::cin >> num_curves1;  
  while (num_curves1--) {
    std::cin >> x1 >> y1 >> x2 >> y2;

    pmwx1.insert(Curve(Point(x1,y1), Point(x2,y2)));
  } 
  
  std::cin >> num_curves2;
  while (num_curves2--) {
    std::cin >> x1 >> y1 >> x2 >> y2;

    pmwx2.insert(Curve(Point(x1,y1), Point(x2,y2))); 
  }
  
  MapOverlay_incremental                 ovl_incremental;

  cout<<"calling map overlay";
  MapOverlay map1(pmwx1, &ovl_incremental), map2(pmwx2, &ovl_incremental);
  MapOverlay map_overlay(map1, map2, &ovl_incremental);

  CGAL::Window_stream W(700, 700);
  W.init(-10, 10, -10);
  W.set_node_width(3);
  
  W.display();
  
  W<<CGAL::BLUE;
  W<<map_overlay.subdivision();
  
  return 0;
}
