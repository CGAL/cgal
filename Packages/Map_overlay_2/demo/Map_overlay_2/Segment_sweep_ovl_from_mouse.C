// demo/Map_overlay_2/Segment_bops_from_mouse.C
//
// constructs an overlay of two segment planar maps from CGAL window.
// We use the leda traits (therefore we are using leda functions).

//#define CGAL_NOTF_DEBUG

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

#include <vector>

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
typedef CGAL::Map_overlay_2<PM, Ovl_change_notification>      MapOverlay;

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

/*void  draw_and_locate_maps (const MapOverlay& ovl, 
                            const PM& pm1, 
                            const PM& pm2, 
                            CGAL::Window_stream& W)
{ 
  
  // debug!
  //cout<<"unbounded faces: "<< &*(ovl.first_creator()->subdivision().unbounded_face())<<
  //  "  "<<&*(ovl.second_creator()->subdivision().unbounded_face())<<endl;

  //if (!(ovl.second_creator()->subdivision().unbounded_face()->is_unbounded()))
  //  cout<<"error: unbounded_face of second map is bounded"<<endl;
  
  W << CGAL::BLUE;
  W << pm1;
  W << CGAL::RED;
  W << pm2;

  for (;;) {
    double  x,y;

    int b = W.read_mouse(x,y);
    if (b == 4) 
      break;

    W.clear();
    W << CGAL::BLUE;
    W << pm1;
    W << CGAL::RED;
    W << pm2;
    
    // debug!
    //MapOverlay map_overlay(pm1, pm2);
    // end debug.
    const PM& pm = ovl.subdivision();
   
    
    NT z(x), w(y);
    Point p(z, w);
    
    PM::Locate_type lt;
    PM::Halfedge_const_handle e = pm.locate(p, lt);
    
    PM::Face_const_handle fh = e->face();
    if (lt == PM::UNBOUNDED_FACE) {
      std::cout << "UNBOUNDED" << endl;
    }
    else {
      W.set_node_width(4);
      
      PM::Face_const_handle f1 = ovl.change_notification()->get_first_face_above(fh);
      PM::Face_const_handle f2 = ovl.change_notification()->get_second_face_above(fh);
      
      //cout<<"f1 and f2: "<<&*f1<<" "<<&*f2<<endl;
      
      //if (f1 == fh || f2 == fh)
      //  std::cout<<"NULL faces pointers\n";
 
      leda_color fg_col;
      if (f1 != fh && !(f1->is_unbounded()) && f2 != fh && !(f2->is_unbounded()) )
        fg_col = leda_violet;
      else if (f1 != fh && !(f1->is_unbounded()))
        fg_col = leda_blue;
      else if (f2 != fh && !(f2->is_unbounded()))
        fg_col = leda_red;
      else
        fg_col = leda_orange;
      
      W.set_color(fg_col);
      W.set_fill_color(fg_col);

      PM::Ccb_halfedge_const_circulator cc = fh->outer_ccb();
      
      leda_list<leda_point> points_list;
      do {
        leda_point p = leda_point(CGAL::to_double(cc->source()->point().xcoordD()), 
                                  CGAL::to_double(cc->source()->point().ycoordD()) );

        points_list.push_back(p);
      } while (++cc != fh->outer_ccb());
      W.draw_filled_polygon(points_list);

      do {
        W << (*cc).curve();
        W << cc->source()->point();
      } while (++cc != fh->outer_ccb());
      
      for (PM::Holes_const_iterator hit = fh->holes_begin(); hit != fh->holes_end(); ++hit) {
        PM::Ccb_halfedge_const_circulator cc(*hit);
        do{
          W << cc->curve();
        } while (++cc != *hit);
      }
      
      for (PM::Vertex_const_iterator  v_iter = pm.vertices_begin(); 
           v_iter != pm.vertices_end(); ++v_iter){
        if (v_iter->get_first_vertex_above() != NULL && v_iter->get_second_vertex_above() != NULL)
          W.set_color(leda_violet);
        else if (v_iter->get_first_vertex_above() != NULL && v_iter->get_second_vertex_above() == NULL)
          W<<CGAL::BLUE;
        else if (v_iter->get_second_vertex_above() != NULL && v_iter->get_first_vertex_above() == NULL)
          W << CGAL::RED;
        else
          W << CGAL::ORANGE;

        W << v_iter->point();
      }

      if (lt == PM::EDGE){
        if (e.operator->()->get_first_halfedge_above() != NULL)
          std::cout<<"edge belongs to first subdivision\n";
        if (e.operator->()->get_second_halfedge_above() != NULL)
          std::cout<<"edge belongs to second subdivision\n";
      }
    }
  }
}*/

int  read_planar_map(PM& pm, CGAL::Window_stream& W)
{
  std::vector<Point> cv1;
  std::vector<Curve> curves;
  Point pnt;
  bool begin=true;
  
  while (true)
    {
      double x, y;
      int b = W.get_mouse(x,y);
      
      if (b==2 || b==3){
        //creating planar map.
        Traits traits;
        CGAL::sweep_to_construct_planar_map_2(curves.begin(), 
                                              curves.end(),
                                              traits, pm);
        W << pm;
        return b;
      }
    
      pnt = Point(x,y);
      
      if (b == MOUSE_BUTTON(1))
        {
          
          for(std::vector<Curve>::iterator iter = curves.begin();
              iter != curves.end(); ++iter) {
            //we are using the leda sqr_dist func
            if ( pnt.sqr_dist(iter->source()) < ((W.xmax() - W.xmin())/50)*((W.xmax() - W.xmin())/50) )
              pnt=iter->source();
            if ( pnt.sqr_dist(iter->target()) < ((W.xmax() - W.xmin())/50)*((W.xmax() - W.xmin())/50) )
              pnt=iter->target();
          }
          
          cv1.push_back(pnt);
          W << CGAL::BLACK;
          W << pnt;
          W << CGAL::GREEN;
          
          if (!begin) {
            if ( cv1[0] == cv1[1] ){
              //Error. Segment has a zero length.
              W.set_status_string("Error. Segment has a zero length.");
              redraw( &W );
            }
            else{ 
              curves.push_back(Curve(cv1[0],cv1[1]));
              W.set_status_string("Left mouse button - segment input. "
                                  "Finish button -query mode");
              W << Curve(cv1[0],cv1[1]);
            }
            cv1.clear();
          }
          begin=!begin;
        }
    }
}

int main()
{
  double x0=-700,x1=700,y0=-700;
  
  cout<<"bogi"<<endl;
  
  W.init(x0,x1,y0);
  W.set_redraw(redraw);
  W.set_mode(leda_src_mode);
  W.set_node_width(3);
  W.open_status_window();
  W.button("Insert Map I",1);
  W.button("Insert Map II",2);
  W.button("Create Overlay",3);
  W.button("Exit",4);
  
  W.display();
  
  //read input from window
  std::cout << "Left button to start and end the segment.\n";
  std::cout << "Clicking close to a vertex assumes the location" 
	    << "is at the vertex"
	    << std::endl;
 
  W.set_status_string( "Left mouse button - segment input. "
		       "Finish button - query mode" );

  while (true){
    double x, y;
    int b = W.get_mouse(x,y);
    if (b==1){
      std::cout<<"Insert first map"<<std::endl;
      b = read_planar_map(pm1,W);
      if (b==2){
        std::cout<<"Insert second map"<<std::endl;
        b=read_planar_map(pm2,W);
      }
    }
    if (b==3)
      break;
  }
  
  W << CGAL::RED;
  W << pm1;
  W << CGAL::BLUE;
  W << pm2;

  //MapOverlay_sweep   ovl_sweep, pm_sweep1, pm_sweep2;
  MapOverlay map1(pm1);
  MapOverlay map2(pm2);
  MapOverlay map_overlay(map1, map2);
  
  //MapOverlay map_overlay(pm1, pm2);  // makes problem with the pointer to the creators: the two contrsucted overlays of pm1 and pm2 are temporary variables, and hence the pointers to them are not valid.
  
  std::cout<<"Locate Overlay Face:"<<endl;
  std::cout<<"Purple Face - an overlay face laying under two bounded faces"<<std::endl;
  std::cout<<"Blue Face - an overlay face laying under a bounded face of the first map and the unbounded face of the second"<<std::endl;
  std::cout<<"Red Face - an overlay face laying under the unbounded face of the first map and a bounded face of the second"<<std::endl;
  std::cout<<"Orange Face - an overlay face laying under the unbounded faces of both maps"<<endl;
  
  // Point Location Queries
  W.set_status_string("Map Overlay.");
  
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
      pm2.halfedges_begin() != pm2.halfedges_end() ){
    CGAL::Ovl_utility<PM,NT> utility;
    utility.draw_and_locate_maps(map_overlay,pm1,pm2,W);
  }
 
  return 0;  
}

#endif // CGAL_USE_LEDA
