// demo/Map_overlay_2/Segment_bops_from_mouse.C
//
// constructs an overlay of two segment planar maps from CGAL window.
// We use the leda traits (therefore we are using leda functions).

#include <CGAL/config.h> // needed for the LONGNAME flag

#ifdef CGAL_CFG_NO_LONGNAME_PROBLEM
// Define shorter names to please linker (g++/egcs)
#define Planar_map_with_intersections Pmwx
#define Arr_leda_segment_exact_traits Alset
#define Bop_default_dcel Bdd
#define Boolean_operation Bo
#define In_place_list_iterator IPLI
#define Arr_2_vertex_base Avb
#define Arr_2_halfedge_base Ahb
#define Arr_2_face_base Afb
#define Point_2 pT
#define Segment_2 sT
#define Topological_map TpM
#define _List_iterator Lit
#define Halfedge hE
#define Forward_circulator_tag Fct
#endif

#include <CGAL/basic.h>
#include <CGAL/Cartesian.h>

#ifndef CGAL_USE_LEDA
int main()
{

  std::cout << "Sorry, this demo needs LEDA for visualisation.";
  std::cout << std::endl;

  return 0;
}

#else

#include <CGAL/Arr_polyline_traits.h>
#include <CGAL/Pm_walk_along_line_point_location.h>
#include <CGAL/Pm_with_intersections.h>


#include <CGAL/Map_overlay_default_dcel.h>
#include <CGAL/Map_overlay.h>
#include <CGAL/Map_overlay_incremental.h>

#include <CGAL/leda_real.h>
#include <LEDA/rat_window.h>
#include <CGAL/IO/Pm_Window_stream.h>
#include <CGAL/Ovl_utility.h>

#if defined(LEDA_NAMESPACE)
using namespace leda;
#endif

typedef leda_real                            NT;
typedef CGAL::Cartesian<NT>                  Rep;
typedef CGAL::Arr_polyline_traits<Rep>       Traits;

typedef Traits::Point_2                      Point;
typedef Traits::X_curve_2                    X_curve;
typedef Traits::Curve_2                      Curve;

typedef CGAL::Map_overlay_default_dcel<Traits>        Dcel;
typedef CGAL::Planar_map_2<Dcel, Traits>              PM;
typedef CGAL::Planar_map_with_intersections_2<PM>     Pmwx;

typedef CGAL::Map_overlay_default_notifier<PM>      MapOverlay_change_notification;
typedef CGAL::Map_overlay_incremental<Pmwx, MapOverlay_change_notification>   
                                                                 MapOverlay_incremental;
typedef CGAL::Map_overlay<Pmwx, MapOverlay_change_notification>  MapOverlay;

typedef CGAL::Pm_walk_along_line_point_location<PM>             PmWalkPL;



// global variables are used so that the redraw function for the LEDA window
// can be defined to draw information found in these variables.
//static PmWalkPL pm_walk1, pm_walk2;
static Pmwx pmwx1; 
static Pmwx pmwx2;
static CGAL::Window_stream W(500, 500, "CGAL - Segment Map-Overlay Demo: Incremental Algorithm");

leda_point to_leda_pnt(Point p)
{
  return leda_point(p.x().to_double(), p.y().to_double());
}

// drawing functions
CGAL::Window_stream& operator<<(CGAL::Window_stream& os,
                          const  Point& p)                         
{
  // conversion to leda_point in order to show it on screen
  return os << to_leda_pnt(p);
}

// draw a polyline, with points as 'x's
CGAL::Window_stream& operator<<(CGAL::Window_stream& os,
				const X_curve &c)
{
  X_curve::const_iterator sit, tit;
  sit = c.begin();
  tit = sit; 
  ++tit;
  
  W.draw_point(to_leda_pnt(*sit), leda_green); // draw first point
  for (; tit != c.end(); ++tit, ++sit)
    {
      // conversion to screen drawble segment
      os << leda_segment(to_leda_pnt(*sit), to_leda_pnt(*tit));
      W.draw_point(to_leda_pnt(*tit), leda_green);
    }
  
  return os;
}

// redraw function for the LEDA window. 
// used automatically when window reappears.
void redraw(CGAL::Window_stream * wp) 
{ wp->start_buffering();
  wp->clear();
  // draw planar map.
  *wp << CGAL::BLUE;
  *wp << pmwx1;
  *wp << CGAL::RED;
  *wp << pmwx2;
  wp->flush_buffer();
  wp->stop_buffering();
}

int main(int argc, char* argv[])
{
  double x0=-200,x1=200,y0=-200;
  
  if (argc != 3) {
    std::cout << "usage: Segment_sweep_ovl_from_file filename\n";
    exit(1);
  }
  
  CGAL::Ovl_utility<Pmwx,NT> utility;
  
  utility.scan_polyline_pmwx(argv[1],pmwx1);
  utility.scan_polyline_pmwx(argv[2],pmwx2);
  
  utility.calc_window_size(pmwx1,x0,x1,y0);
  utility.calc_window_size(pmwx2,x0,x1,y0);
  
  W.init(x0,x1,y0);
  W.set_redraw(redraw);
  W.set_mode(leda_src_mode);
  W.set_node_width(3);
  W.open_status_window();
  W.button("Exit",4);
  W.display();
  
  //POINT LOCATION
  W.set_status_string( "Enter a point with left button."); 
  
  W << CGAL::RED;
  W << pmwx1;
  W << CGAL::BLUE;
  W << pmwx2;
  
  MapOverlay_incremental   ovl_incremental, pmwx_incremental1, pmwx_incremental2;
  MapOverlay map1(pmwx1, &pmwx_incremental1);
  MapOverlay map2(pmwx2, &pmwx_incremental2);
  MapOverlay map_overlay(map1, map2, &ovl_incremental);
  
  // Point Location Queries
  W.set_status_string("Map Overlay. ");
  
  std::cout<<"Locate Overlay Face:"<<endl;
  std::cout<<"Purple Face - an overlay face laying under two bounded faces"<<std::endl;
  std::cout<<"Blue Face - an overlay face laying under a bounded face of the first map and the unbounded face of the second"<<std::endl;
  std::cout<<"Red Face - an overlay face laying under the unbounded face of the first map and a bounded face of the second"<<std::endl;
  std::cout<<"Orange Face - an overlay face laying under the unbounded faces of both maps"<<endl;

  // if first map is empty
  if (pmwx1.halfedges_begin() == pmwx1.halfedges_end()) 
    {
      std::cout << std::endl;
      std::cout << "No edges were inserted to the first planar map. First Planar map is empty. Exiting.";
      std::cout << std::endl;
    }
  
  // if second map is empty
  if (pmwx2.halfedges_begin() == pmwx2.halfedges_end()) 
    {
      std::cout << std::endl;
      std::cout << "No edges were inserted to the first planar map. First Planar map is empty. Exiting.";
      std::cout << std::endl;
    }
  
  if (pmwx1.halfedges_begin() != pmwx1.halfedges_end() && 
      pmwx2.halfedges_begin() != pmwx2.halfedges_end() )
    utility.draw_and_locate_maps(map_overlay,pmwx1,pmwx2,W);
 
  return 0;  
}

#endif // CGAL_USE_LEDA
