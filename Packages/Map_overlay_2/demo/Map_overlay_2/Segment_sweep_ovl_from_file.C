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

#include <LEDA/rat_window.h>
#include <CGAL/IO/Pm_Window_stream.h>
#include <vector>

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
static PmWalkPL pm_walk1, pm_walk2;
static PM pm1(&pm_walk1); 
static PM pm2(&pm_walk2);
static CGAL::Window_stream W(700, 700, "CGAL - Segment Map-Overlay Demo: Sweep Algorithm");

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
  W << CGAL::BLUE;
  *wp << pm1; 
  W << CGAL::RED;
  *wp << pm2;
  
  wp->flush_buffer();
  wp->stop_buffering();
}

void  draw_and_locate_maps (const MapOverlay& ovl, 
                            const PM& pm1, 
                            const PM& pm2, 
                            CGAL::Window_stream& W)
{ 
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
      
      if (f1 == fh || f2 == fh)
        std::cout<<"NULL faces pointers\n";
 
      leda_color fg_col;
      if (f1 != fh && !(f1->is_unbounded()) && f2 != fh && !(f2->is_unbounded()) ){
        // debug!
        cout<<"both maps lay above face"<<endl;
        fg_col = leda_violet;
      }
      else if (f1 != fh && !(f1->is_unbounded())){
        cout<<"first map lays above face"<<endl;
        fg_col = leda_blue;
      }
      else if (f2 != fh && !(f2->is_unbounded())){
        cout<<"second map lays above face"<<endl;
        fg_col = leda_red;
      }
      else{
        fg_col = leda_orange;
        cout<<"non map lays above face"<<endl;
      }
      
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
        std::cout<<cc->curve()<<std::endl;
      } while (++cc != fh->outer_ccb());
      
      for (PM::Holes_const_iterator hit = fh->holes_begin(); hit != fh->holes_end(); ++hit, ++hit) {
        PM::Ccb_halfedge_const_circulator cc(*hit);
        do{
        W << cc->curve();
        } while (++cc != *hit);
      }
      
      for (PM::Vertex_const_iterator  v_iter = pm.vertices_begin(); v_iter !=  pm.vertices_end(); v_iter++){
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
}

void  scan_planar_map(const char* filename, PM& pm)
{
  std::vector<Curve> curves;
  
  std::ifstream file(filename);
  int num_curves;
  file >> num_curves;
  while (num_curves--) {
    //NT x,y;
    double x,y; //(actually int)
    file >> x >> y;
    NT nx(x),ny(y);
    Point s(nx,ny);
    file >> x >> y;
    NT mx(x),my(y);
    Point t(mx,my);
    
    curves.push_back(Curve(s,t));
  }
  
  Traits traits;
  CGAL::sweep_to_construct_planar_map_2(curves.begin(), 
                                        curves.end(),
                                        traits, pm);
  W << pm;
  
}
  
int main(int argc, char* argv[])
{
  double x0=-200,x1=200,y0=-200;
  
  if (argc != 3) {
    std::cout << "usage: Segment_sweep_ovl_from_file filename\n";
    exit(1);
  }
  
  W.init(x0,x1,y0);
  W.set_redraw(redraw);
  W.set_mode(leda_src_mode);
  W.set_node_width(3);
  W.open_status_window();
  W.button("Exit",4);
  
  W.display();
  
  scan_planar_map(argv[1],pm1);
  scan_planar_map(argv[2],pm2);
 
  //POINT LOCATION
  W.set_status_string( "Enter a point with left button.");
  
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
  W.set_status_string("Map Overlay. ");
  
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
    draw_and_locate_maps(map_overlay,pm1,pm2,W);
 
  return 0;  
}

#endif // CGAL_USE_LEDA


