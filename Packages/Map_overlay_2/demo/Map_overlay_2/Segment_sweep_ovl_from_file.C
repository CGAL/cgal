// demo/Map_overlay_2/Segment_bops_from_mouse.C
//
// constructs an overlay of two segment planar maps from CGAL window.
// We use the leda traits (therefore we are using leda functions).

//#define CGAL_NOTF_DEBUG

#include <CGAL/config.h> // needed for the LONGNAME flag

#ifdef CGAL_CFG_NO_LONGNAME_PROBLEM
// Define shorter names to please linker (g++/egcs)
#define Planar_map Pm
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

#include <CGAL/sweep_to_construct_planar_map_2.h>

#ifndef CGAL_MAP_OVERLAY_DEFAULT_DCEL_H
#include <CGAL/Map_overlay_default_dcel.h>
#endif

#ifndef CGAL_MAP_OVERLAY_H
#include <CGAL/Map_overlay.h>
#endif

#include <CGAL/leda_rational.h>
#include <LEDA/rat_window.h>
#include <CGAL/IO/Pm_Window_stream.h>
#include <CGAL/Ovl_utility.h>

#if defined(LEDA_NAMESPACE)
using namespace leda;
#endif

typedef leda_rational                               NT;
typedef CGAL::Arr_leda_segment_exact_traits         Traits;

typedef Traits::Point_2                             Point;
typedef Traits::Curve_2                             Curve;
typedef Traits::X_curve                             X_curve;

typedef CGAL::Map_overlay_default_dcel<Traits>      Dcel;
typedef CGAL::Planar_map_2<Dcel, Traits>            PM;

typedef CGAL::Map_overlay_default_notifier<PM>                Ovl_change_notification;
typedef CGAL::Map_overlay<PM, Ovl_change_notification>        MapOverlay;

typedef CGAL::Pm_walk_along_line_point_location<PM>           PmWalkPL;

// global variables are used so that the redraw function for the LEDA window
// can be defined to draw information found in these variables.
static PmWalkPL pm_walk1, pm_walk2;
static PM pm1(&pm_walk1); 
static PM pm2(&pm_walk2);
static CGAL::Window_stream W(500, 500, "CGAL - Segment Map-Overlay Demo: Sweep Algorithm");


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
  double x0=-200,x1=200,y0=-200;
  
  if (argc != 3) {
    std::cout << "usage: Segment_sweep_ovl_from_file filename\n";
    exit(1);
  }
  
  CGAL::Ovl_utility<PM,NT> utility;
  
  utility.scan_segment_planar_map(argv[1],pm1);
  utility.scan_segment_planar_map(argv[2],pm2);
 
  utility.calc_window_size(pm1,x0,x1,y0);
  utility.calc_window_size(pm2,x0,x1,y0);
  
  W.init(x0,x1,y0);
  W.set_redraw(redraw);
  W.set_mode(leda_src_mode);
  W.set_node_width(3);
  W.open_status_window();
  W.button("Exit",4);
  W.display();
  
  W << CGAL::RED;
  W << pm1;
  W << CGAL::BLUE;
  W << pm2;

  //MapOverlay_sweep   ovl_sweep, pm_sweep1, pm_sweep2;
  PmWalkPL ovl_walk;
  MapOverlay map1(pm1);
  MapOverlay map2(pm2);
  MapOverlay map_overlay(map1, map2, &ovl_walk);
  
  //MapOverlay map_overlay(pm1, pm2);  // makes problem with the pointer to the creators: the two contrsucted overlays of pm1 and pm2 are temporary variables, and hence the pointers to them are not valid.
  
  std::cout<<"Locate Overlay Face:"<<endl;
  std::cout<<"Purple Face - an overlay face laying under two bounded faces"<<std::endl;
  std::cout<<"Blue Face - an overlay face laying under a bounded face of the first map and the unbounded face of the second"<<std::endl;
  std::cout<<"Red Face - an overlay face laying under the unbounded face of the first map and a bounded face of the second"<<std::endl;
  std::cout<<"Orange Face - an overlay face laying under the unbounded faces of both maps"<<endl;
  
  // Point Location Queries
  W.set_status_string("Map Overlay. Enter a point with left button.");
  
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
    utility.draw_and_locate_maps(map_overlay,pm1,pm2,W);
 
  return 0;  
}

#endif // CGAL_USE_LEDA
