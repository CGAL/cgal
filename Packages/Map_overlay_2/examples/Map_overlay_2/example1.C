#include <CGAL/Cartesian.h>  //CGAL definitions that need to be before anything else
#include <CGAL/Quotient.h>
#include <CGAL/Pm_walk_along_line_point_location.h>
#include <CGAL/Map_overlay_default_dcel.h>
#include <CGAL/Map_overlay.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Planar_map_2.h>
#include <CGAL/IO/Pm_iostream.h>

#include <iostream>
#include <vector>
#include <list>

// uncomment if LEDA is installed.
//#include <CGAL/IO/cgal_window.h>  //used for visualization.
//#include <CGAL/IO/Pm_Window_stream.h>

typedef CGAL::Quotient<int>                                     NT;
typedef CGAL::Cartesian<NT>                                     K;
typedef CGAL::Arr_segment_traits_2<K>                           Traits;
typedef Traits::Point_2                                         Point_2;
typedef Traits::X_curve_2                                       X_curve_2;
typedef Traits::Curve_2                                         Curve_2;

typedef CGAL::Map_overlay_default_dcel<Traits>                  Dcel;
typedef CGAL::Planar_map_2<Dcel,Traits>                         Planar_map;
typedef CGAL::Planar_map_with_intersections_2<Planar_map>       Pmwx;
typedef CGAL::Map_overlay_2<Planar_map>                         MapOverlay;

typedef CGAL::Pm_walk_along_line_point_location<Planar_map>     PmWalkPL;

using std::cin;
using std::cout;
using std::endl;


int main()
{
  std::vector<Curve_2> segments;
  PmWalkPL     pl_walk1, pl_walk2;
  Pmwx pm1(&pl_walk1), pm2(&pl_walk2);
  int   num_curves1, num_curves2;
  
  std::cin >> num_curves1;
  
  NT x1, y1, x2, y2;
  
  while (num_curves1--) {
    std::cin >> x1 >> y1 >> x2 >> y2;

    Point_2 p1(x1,y1), p2(x2,y2);
    segments.push_back(Curve_2(p1, p2));
  } 
  
  pm1.insert(segments.begin(), segments.end());

  segments.clear();
  
  std::cin >> num_curves2;
  
  while (num_curves2--) {
    std::cin >> x1 >> y1 >> x2 >> y2;
    
    Point_2 p1(x1,y1), p2(x2,y2);
    segments.push_back(Curve_2(p1, p2)); 
  }
  
  pm2.insert(segments.begin(), segments.end());
  
  cout << "Calling map overlay" << endl;
  cout << endl;
  
  PmWalkPL  ovl_walk;
  //MapOverlay map1(pm1), map2(pm2);
  MapOverlay map_overlay(pm1, pm2, &ovl_walk);

  cout << "Writing the resulting subdivision induced the the overlay of the two input subdivision" << endl;
  cout << endl;

  std::cout<<map_overlay.subdivision();
  
  // uncomment if LEDA is installed.
  //CGAL::Window_stream W(700, 700);
  //W.init(-10, 10, -10);
  //W.set_node_width(3);
  //W.display();
  //W<<CGAL::BLUE;
  //W<<map_overlay.subdivision();
  
  return 0;
}


