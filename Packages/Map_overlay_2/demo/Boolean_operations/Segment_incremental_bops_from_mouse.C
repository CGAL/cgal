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

#include <vector>

#include <CGAL/Arr_leda_segment_exact_traits.h>
#include <CGAL/Pm_walk_along_line_point_location.h>
#include <CGAL/Pm_with_intersections.h>

#include <CGAL/Bop_default_dcel.h>
#include <CGAL/Map_overlay.h>
#include <CGAL/Map_overlay_incremental.h>
#include <CGAL/Boolean_operations_2.h>

#include <CGAL/leda_rational.h>
#include <LEDA/rat_window.h>
#include <CGAL/IO/Pm_Window_stream.h>
#include <CGAL/Bops_utility.h>


#if defined(LEDA_NAMESPACE)
using namespace leda;
#endif

typedef leda_rational                                NT;
typedef CGAL::Arr_leda_segment_exact_traits          Traits;

typedef Traits::Point                                Point;
typedef Traits::X_curve                              X_curve;

typedef CGAL::Bop_default_dcel<Traits>               Dcel;
typedef CGAL::Planar_map_2<Dcel, Traits>             PM;
typedef CGAL::Planar_map_with_intersections_2<PM>    Pmwx;

typedef CGAL::Map_overlay_default_notifier<PM>      MapOverlay_change_notification;
typedef CGAL::Map_overlay_incremental<Pmwx, 
              MapOverlay_change_notification>   MapOverlay_incremental;
typedef CGAL::Map_overlay_2<Pmwx, MapOverlay_change_notification>  
                                                          MapOverlay;
typedef CGAL::Boolean_operations_2<MapOverlay>                  Bops;

typedef CGAL::Pm_walk_along_line_point_location<PM>       PmWalkPL;

typedef Bops::Faces_container                      Faces_container;
typedef Bops::Halfedges_container                  Halfedges_container;
typedef Bops::Vertices_container                   Vertices_container;

// global variables are used so that the redraw function for the LEDA window
// can be defined to draw information found in these variables.
//static PmWalkPL pm_walk1, pm_walk2;
static Pmwx pmwx1; 
static Pmwx pmwx2;
static CGAL::Window_stream W(700, 700, "CGAL - Segment Boolean-Operations Demo");

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

/*void  draw_and_locate_maps (Bops& bops , 
                            const Pmwx& arr1, 
                            const Pmwx& arr2, 
                            CGAL::Window_stream& W)
{
  MapOverlay_change_notification   tmp_notf;
  
  W<<CGAL::BLUE;
  W<<arr1;
  W<<CGAL::RED;
  W<<arr2;

  for (;;) {
    double  x,y;

    int b = W.read_mouse(x,y);
    if (b == 7) 
      break;
    
    W.set_node_width(3);
    W.clear();
    W<<CGAL::BLUE;
    W<<arr1;
    W<<CGAL::RED;
    W<<arr2;

    leda_color fg_col = leda_violet;
    W.set_color(fg_col);
    W.set_fill_color(fg_col);

    if (b == 4 || b == 5 || b == 6){
      Faces_container       faces_result;
      Halfedges_container   halfedges_result;
      Vertices_container    vertices_result;

      if (b == 4)
        bops.intersection(faces_result,halfedges_result,vertices_result);

      else if (b == 5)
        bops.Union(faces_result,halfedges_result,vertices_result);
      
      else if (b == 6)
        bops.symmetric_difference(faces_result, halfedges_result, vertices_result);

      std::cout<<"number of faces is "<< faces_result.size()<<std::endl;
      std::cout<<"The number of resulting halfedges is "<<halfedges_result.size()<<std::endl;
      std::cout<<"The number of resulting vertices is "<<vertices_result.size()<<std::endl;

      for (list<Pmwx::Face_const_handle>::iterator  f_iter = faces_result.begin(); 
           f_iter != faces_result.end(); f_iter++){
        Pmwx::Face_const_handle fh = *f_iter;
        if (fh->is_unbounded())
          continue;
        
        W.set_color(fg_col);
        W.set_fill_color(fg_col);
        
        Pmwx::Ccb_halfedge_const_circulator cc = fh->outer_ccb();
        leda_list<leda_point> points_list;
        do {
          leda_point p = cc->source()->point().to_point(); 
          
          points_list.push_back(p);
        } while (++cc != fh->outer_ccb());
        W.draw_filled_polygon(points_list);
        
        for (Pmwx::Holes_const_iterator hit = fh->holes_begin(); hit != fh->holes_end(); ++hit) {
          leda_list<leda_point> holes_points_list;
          Pmwx::Ccb_halfedge_const_circulator cc(*hit);
          W.set_fill_color(leda_black);
          W << CGAL::WHITE;
          do{
            W << cc->curve();
            
            leda_point p = cc->source()->point().to_point(); 
            points_list.push_back(p);
          } while (++cc != *hit);
          W.draw_filled_polygon(holes_points_list);
        }
        
        do {
          W <<CGAL::WHITE;
          W << (*cc).curve();
          //std::cout<<cc->curve()<<std::endl;
        } while (++cc != fh->outer_ccb());
      }
      
      for (list<Pmwx::Halfedge_const_handle>::iterator h_iter = halfedges_result.begin(); 
           h_iter != halfedges_result.end(); h_iter++){
        Pmwx::Halfedge_const_handle h = *h_iter;

        W.set_color(fg_col);
        W.set_node_width(4);
        W << h->curve();
      }
      
      for (list<Pmwx::Vertex_const_handle>::iterator  v_iter = vertices_result.begin(); 
           v_iter != vertices_result.end(); v_iter++){
        W.set_node_width(4);
        //for (Vertex_const_iterator  v_iter = bop.get_map_overlay().get_arr().vertices_begin(); 
        //v_iter !=  bop.get_map_overlay().get_arr().vertices_end(); v_iter++){
        
        Pmwx::Vertex_const_handle vh = *v_iter;
        //Vertex_const_iterator vh = v_iter;
        
        W.set_color(leda_violet);
        W << vh->point();
      }
    }
  }
}*/


int  read_pmwx(Pmwx& pmwx, CGAL::Window_stream& W)
{
  std::vector<Point> cv1;
  Point pnt;
  bool begin=true;
  
  while (true)
    {
      double x, y;
      int b = W.get_mouse(x,y);
      
      if (b==2 || b ==3) 
        return b;
    
      pnt = Point(x,y);
      
      if (b == MOUSE_BUTTON(1))
        {
          
          for(Pmwx::Vertex_iterator vi = pmwx.vertices_begin();
              vi != pmwx.vertices_end(); ++vi) {
            //we are using the leda sqr_dist func
            if ( pnt.sqr_dist(vi->point()) < ((W.xmax() - W.xmin())/50)*((W.xmax() - W.xmin())/50) )
              pnt=vi->point();
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
              pmwx.insert(X_curve(cv1[0],cv1[1]));
              W.set_status_string("Left mouse button - segment input. "
                                  "Finish button -query mode");
              W << pmwx;
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
  
  W.init(x0,x1,y0);
  W.set_redraw(redraw);
  W.set_mode(leda_src_mode);
  W.set_node_width(3);
  W.open_status_window();
  W.button("Insert Map I",1);
  W.button("Insert Map II",2);
  W.button("Start Bops",3);
  W.button("Intersection",4);
  W.button("Union",5);
  W.button("Symmetric Difference",6);
  W.button("Exit",7);
  W.display();
  
  //read input from window
  std::cout << "Left button to start and end the segment.\n";
  std::cout << "Clicking close to a vertex assumes the location" 
	    << "is at the vertex"
	    << std::endl;
 
  W.set_status_string( "Left mouse button - segment input. "
		       "Finish button -query mode" );
  
  
  while (true){
     double x, y;
     int b = W.get_mouse(x,y);
     if (b==1){
       std::cout<<"Insert first map"<<std::endl;
       b = read_pmwx(pmwx1,W);
       if (b==2){
         std::cout<<"Insert second map"<<std::endl;
         b=read_pmwx(pmwx2,W);
       }
     }
     if (b==3)
       break;
  }
  
  pmwx1.unbounded_face()->set_ignore_bop(false); 
  pmwx2.unbounded_face()->set_ignore_bop(false);
  
  W << CGAL::RED;
  W << pmwx1;
  W << CGAL::BLUE;
  W << pmwx2;
  
  MapOverlay_incremental     ovl_incremental, pmwx_incremental1, pmwx_incremental2;
  MapOverlay map1(pmwx1, &pmwx_incremental1);
  MapOverlay map2(pmwx2, &pmwx_incremental2);
  Bops bops(MapOverlay(map1, map2, &ovl_incremental));
  
  // Point Location Queries
  W.set_status_string("Boolean Operations. "
		      "Finish button - exit." );
  
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
      pmwx2.halfedges_begin() != pmwx2.halfedges_end() ){
    CGAL::Bops_utility<Pmwx,NT> utility;
    utility.draw_and_locate_maps(bops,pmwx1,pmwx2,W);
  }
 
  return 0;  
}

#endif // CGAL_USE_LEDA
