#include <CGAL/basic.h> //CGAL definitions that need to be before anything else
#include <CGAL/Cartesian.h>
#include <CGAL/Quotient.h>
#include <CGAL/Pm_walk_along_line_point_location.h>
#include <CGAL/Map_overlay_default_dcel.h>
#include <CGAL/Map_overlay.h>
#include <CGAL/Arr_segment_exact_traits.h>
#include <CGAL/Planar_map_2.h>
#include <CGAL/sweep_to_construct_planar_map_2.h>
#include <CGAL/IO/Pm_iostream.h>

#include <iostream.h>
#include <vector>
#include <list>

// uncomment if LEDA is installed.
//#include <CGAL/IO/cgal_window.h>  //used for visualization.
//#include <CGAL/IO/Pm_Window_stream.h>

typedef CGAL::Quotient<int>                 NT;
typedef CGAL::Cartesian<NT>                 R;
typedef CGAL::Arr_segment_exact_traits<R>   Traits;
typedef Traits::Point                       Point;
typedef Traits::X_curve                     X_curve;
typedef Traits::Curve                       Curve;

typedef CGAL::Map_overlay_default_dcel<Traits> Dcel;
typedef CGAL::Planar_map_2<Dcel,Traits>        Planar_map;
typedef CGAL::Map_overlay_2<Planar_map>        MapOverlay;

typedef CGAL::Pm_walk_along_line_point_location<Planar_map>  PmWalkPL;

using std::cin;
using std::cout;
using std::endl;


int  main()
{
  std::vector<Curve> segments;
  PmWalkPL     pl_walk1, pl_walk2;
  Planar_map   pm1(&pl_walk1), pm2(&pl_walk2);
  int   num_curves1, num_curves2;
  
  std::cin >> num_curves1;
  
  NT          x1, y1, x2, y2;
  
  while (num_curves1--) {
    std::cin >> x1 >> y1 >> x2 >> y2;
    segments.push_back(Curve(Point(x1,y1), Point(x2,y2)));
  } 
  
  Traits traits;
  CGAL::sweep_to_construct_planar_map_2(segments.begin(),segments.end(), 
                                        traits, pm1);

  segments.clear();
  
  std::cin >> num_curves2;
  
  while (num_curves2--) {
    std::cin >> x1 >> y1 >> x2 >> y2;

    segments.push_back(Curve(Point(x1,y1), Point(x2,y2))); 
  }
  
  CGAL::sweep_to_construct_planar_map_2(segments.begin(),segments.end(), 
                                        traits, pm2);
  
  cout<<"Calling map overlay"<<endl;
  cout<<endl;
  
  PmWalkPL  ovl_walk;
  MapOverlay map1(pm1), map2(pm2);
  MapOverlay map_overlay(map1, map2, &ovl_walk);

  cout<<"Writing the resulting subdivision induced the the overlay of the two input subdivision"<<endl;
  cout<<endl;

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
