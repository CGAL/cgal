// demo/Map_overlay_2/Segment_bops_from_mouse.C
//
// constructs an overlay of two segment planar maps from CGAL window.
// We use the leda traits (therefore we are using leda functions).

#include "short_names.h"

#include <CGAL/basic.h>

#ifndef CGAL_USE_LEDA
int main()
{

  std::cout << "Sorry, this demo needs LEDA for visualisation.";
  std::cout << std::endl;

  return 0;
}

#else

#include <CGAL/Arr_leda_segment_exact_traits.h>
#include <CGAL/Pm_walk_along_line_point_location.h>
#include <CGAL/Planar_map_2.h>

#include <CGAL/Arr_segment_circle_traits.h>
#include <CGAL/IO/Segment_circle_Window_stream.h>
#include <CGAL/Bop_default_dcel.h>
#include <CGAL/Map_overlay.h>
#include <CGAL/Map_overlay_naive.h>
#include <CGAL/Boolean_operations_2.h>

#include <CGAL/leda_real.h>
#include <LEDA/rat_window.h>
#include <CGAL/IO/Pm_Window_stream.h>
#include <CGAL/Bops_utility.h>

#if defined(LEDA_NAMESPACE)
using namespace leda;
#endif

typedef leda_real                               NT;
typedef CGAL::Arr_segment_circle_traits<NT>     Traits; 

typedef Traits::Point_2                         Point;
typedef Traits::Segment_2                       Segment;
typedef Traits::Circle                          Circle;
typedef Traits::Conic                           Conic;
typedef Traits::Curve_2                         Curve; 
typedef Traits::X_curve_2                       X_curve;

typedef CGAL::Bop_default_dcel<Traits>          Dcel;
typedef CGAL::Planar_map_2<Dcel, Traits>        PM;

typedef CGAL::Map_overlay_default_notifier<PM>  Ovl_change_notification;
typedef CGAL::Map_overlay_2<PM, Ovl_change_notification>        MapOverlay;
typedef CGAL::Boolean_operations_2<MapOverlay>  Bops;

typedef CGAL::Pm_walk_along_line_point_location<PM>             PmWalkPL;

// global variables are used so that the redraw function for the LEDA window
// can be defined to draw information found in these variables.
static PmWalkPL pm_walk1, pm_walk2;
static PM pm1(&pm_walk1); 
static PM pm2(&pm_walk2);
static CGAL::Window_stream W(500, 500, "CGAL - Segment Boolean-Operations Demo");

// redraw function for the LEDA window. 
// used automatically when window reappears.
void redraw(CGAL::Window_stream * wp) 
{ wp->start_buffering();
  wp->clear();
  // draw planar map.
  W << CGAL::BLUE;
  *wp << pm1; 
  W << CGAL::RED;
  *wp << pm2;
  
  wp->flush_buffer();
  wp->stop_buffering();
}

int main(int argc, char* argv[])
{
  double x0=-700,x1=700,y0=-700;
  
  if (argc != 3) {
    std::cout << "usage: Segment_sweep_ovl_from_file filename\n";
    exit(1);
  }
  
  CGAL::Bops_utility<PM,NT> utility;
  
  utility.scan_seg_circ_planar_map(argv[1],pm1);
  utility.scan_seg_circ_planar_map(argv[2],pm2);

  pm1.unbounded_face()->set_ignore_bop(false); 
  pm2.unbounded_face()->set_ignore_bop(false);
  
  PmWalkPL   ovl_walk;
  MapOverlay map1(pm1);
  MapOverlay map2(pm2);
  Bops bops(MapOverlay(map1, map2, &ovl_walk));
  //Bops bops(pm1, pm2);
  
  utility.calc_window_size(pm1, x0, x1, y0);
  utility.calc_window_size(pm2, x0, x1, y0);
  
  W.init(x0,x1,y0);
  W.set_redraw(redraw);
  W.set_mode(leda_src_mode);
  W.set_node_width(3);
  W.open_status_window();
  W.button("Intersection",4);
  W.button("Union",5);
  W.button("Symmetric Difference",6);
  W.button("Exit",7);
  W.display();
   
  W << CGAL::RED;
  W << pm1;
  W << CGAL::BLUE;
  W << pm2;

  // Point Location Queries
  W.set_status_string("Boolean Operations. "
		      "Finish button - Exit." );
  
  // if first map is empty
  if (pm1.halfedges_begin() == pm1.halfedges_end()) 
    {
      std::cout << std::endl;
      std::cout << "No edges were inserted to the first planar map. First Planar map is empty. Exiting.";
      std::cout << std::endl;
    }
  
  // if second map is empty
  if (pm2.halfedges_begin() == pm2.halfedges_end()) 
    {
      std::cout << std::endl;
      std::cout << "No edges were inserted to the first planar map. First Planar map is empty. Exiting.";
      std::cout << std::endl;
    }
  
  if (pm1.halfedges_begin() != pm1.halfedges_end() && 
      pm2.halfedges_begin() != pm2.halfedges_end() )
    utility.draw_and_locate_maps(bops,pm1,pm2,W);
 
  return 0;  
}

#endif  // CGAL_USE_LEDA

