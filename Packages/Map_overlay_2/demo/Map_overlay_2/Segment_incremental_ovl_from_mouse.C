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
#include <CGAL/Pm_with_intersections.h>


#include <CGAL/Map_overlay_default_dcel.h>
#include <CGAL/Map_overlay.h>
#include <CGAL/Map_overlay_incremental.h>


#include <LEDA/rat_window.h>
#include <CGAL/IO/Pm_Window_stream.h>
#include <vector>

#if defined(LEDA_NAMESPACE)
using namespace leda;
#endif

typedef leda_rational                                  NT;
typedef CGAL::Arr_leda_segment_exact_traits            Traits;

typedef Traits::Point                                  Point;
typedef Traits::X_curve                                X_curve;

typedef CGAL::Map_overlay_default_dcel<Traits>                Dcel;
typedef CGAL::Planar_map_2<Dcel, Traits>              PM;
typedef CGAL::Planar_map_with_intersections_2<PM>     Pmwx;

typedef CGAL::Map_overlay_default_notifier<PM>      MapOverlay_change_notification;
typedef CGAL::Map_overlay_incremental<Pmwx, MapOverlay_change_notification>   
                                                                 MapOverlay_incremental;
typedef CGAL::Map_overlay<Pmwx, MapOverlay_change_notification>  MapOverlay;

typedef CGAL::Pm_walk_along_line_point_location<PM>             PmWalkPL;


//I had to add these in global namespace for the program to compile

/*
  CGAL::Window_stream& operator<<(CGAL::Window_stream& os,
  const Point& p)
  {
  //return os << leda_point(p.xcoordD(),p.ycoordD());
  return os << p.to_point();
  }
  
  
  CGAL::Window_stream& operator<<(CGAL::Window_stream& os,
  const X_curve &c)
  {
  return os << c.to_segment();
  }*/

// global variables are used so that the redraw function for the LEDA window
// can be defined to draw information found in these variables.
//static PmWalkPL pm_walk1, pm_walk2;
static Pmwx pmwx1; 
static Pmwx pmwx2;
static CGAL::Window_stream W(700, 700, "CGAL - Segment Arrangement Demo");

/*CGAL_BEGIN_NAMESPACE
  Window_stream& operator<<(Window_stream& os, const PM &pm)
  {
  My_Arr_drawer< PM,
  PM::Ccb_halfedge_const_circulator, 
  PM::Holes_const_iterator> drawer(os);
  
  draw_pm(pm, drawer, os);
  
  return os;
  }
  CGAL_END_NAMESPACE*/

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

void  draw_and_locate_maps (const Pmwx& pmwx, 
                            const Pmwx& pmwx1, 
                            const Pmwx& pmwx2, 
                            CGAL::Window_stream& W)
{ 
  W << CGAL::BLUE;
  W << pmwx1;
  W << CGAL::RED;
  W << pmwx2;

  for (;;) {
    double  x,y;

    int b = W.read_mouse(x,y);
    if (b == 4) 
      break;

    W.clear();
    W << CGAL::BLUE;
    W << pmwx1;
    W << CGAL::RED;
    W << pmwx2;

    //W << pm;
    
    NT z(x), w(y);
    Point p(z, w);

    //debug.
    cout<<p<<endl;
    
    Pmwx::Locate_type lt;
    Pmwx::Halfedge_const_handle e = pmwx.locate(p, lt);
    
    Pmwx::Face_const_handle fh = e->face();
    if (lt == Pmwx::UNBOUNDED_FACE) {
      std::cout << "UNBOUNDED" << endl;
    }
    else {
      W.set_node_width(4);
      
      Pmwx::Face* f1 = (Pmwx::Face*) fh->get_first_face_above();
      Pmwx::Face* f2 = (Pmwx::Face*) fh->get_second_face_above();

      if (f1 == NULL || f2 == NULL)
        std::cout<<"NULL faces pointers\n";
 
      leda_color fg_col;
      if (f1 != NULL && !(f1->is_unbounded()) && f2 != NULL && !(f2->is_unbounded()) ){
        // debug!
        cout<<"both maps lay above face"<<endl;
        fg_col = leda_violet;
      }
      else if (f1 != NULL && !(f1->is_unbounded())){
        cout<<"first map lays above face"<<endl;
        fg_col = leda_blue;
      }
      else if (f2 != NULL && !(f2->is_unbounded())){
        cout<<"second map lays above face"<<endl;
        fg_col = leda_red;
      }
      else{
        fg_col = leda_orange;
        cout<<"non map lays above face"<<endl;
      }
      
      W.set_color(fg_col);
      W.set_fill_color(fg_col);

      Pmwx::Ccb_halfedge_const_circulator cc = fh->outer_ccb();
      
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
        std::cout<<cc->curve()<<std::endl;
      } while (++cc != fh->outer_ccb());
      
      for (Pmwx::Holes_const_iterator hit = fh->holes_begin(); hit != fh->holes_end(); ++hit, ++hit) {
        Pmwx::Ccb_halfedge_const_circulator cc(*hit);
        do{
        W << cc->curve();
        //std::cout<<cc->curve()<<std::endl;
        } while (++cc != *hit);
      }
      
      for (Pmwx::Vertex_const_iterator  v_iter = pmwx.vertices_begin(); v_iter !=  pmwx.vertices_end(); v_iter++){
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

      if (lt == Pmwx::EDGE){
        if (e.operator->()->get_first_halfedge_above() != NULL)
          std::cout<<"edge belongs to first subdivision\n";
        if (e.operator->()->get_second_halfedge_above() != NULL)
          std::cout<<"edge belongs to second subdivision\n";
      }
    }
  }
}

int  read_planar_map(Pmwx& pmwx, CGAL::Window_stream& W)
{
  std::vector<Point> cv1;
  Point pnt;
  bool begin=true;
  
  while (true)
    {
      double x, y;
      int b = W.get_mouse(x,y);
      
      if (b==2 || b==3) 
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
            cout<<"bogi1"<<endl;
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
		       "Finish button -query mode" );
  
  while (true){
    double x, y;
    int b = W.get_mouse(x,y);
    if (b==1){
      std::cout<<"Insert first map"<<std::endl;
      b = read_planar_map(pmwx1,W);
      if (b==2){
        std::cout<<"Insert second map"<<std::endl;
        b=read_planar_map(pmwx2,W);
      }
    }
    if (b==3)
      break;
  }
  
  W << CGAL::RED;
  W << pmwx1;
  W << CGAL::BLUE;
  W << pmwx2;
  
  std::cout<<"Locate Overlay Face:"<<endl;
  std::cout<<"Purple Face - an overlay face laying under two bounded faces"<<std::endl;
  std::cout<<"Blue Face - an overlay face laying under a bounded face of the first map and the unbounded face of the second"<<std::endl;
  std::cout<<"Red Face - an overlay face laying under the unbounded face of the first map and a bounded face of the second"<<std::endl;
  std::cout<<"Orange Face - an overlay face laying under the unbounded faces of both maps"<<endl;
  
  MapOverlay_incremental   ovl_incremental, pmwx_incremental1, pmwx_incremental2;
  MapOverlay map1(pmwx1, &pmwx_incremental1);
  MapOverlay map2(pmwx2, &pmwx_incremental2);
  MapOverlay map_overlay(map1, map2, &ovl_incremental);
  
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
      pmwx2.halfedges_begin() != pmwx2.halfedges_end() )
    draw_and_locate_maps(map_overlay.subdivision(),pmwx1,pmwx2,W);
 
  return 0;  
}

#endif // CGAL_USE_LEDA



