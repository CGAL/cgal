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

#include <CGAL/Bop_default_dcel.h>
#include <CGAL/Map_overlay.h>
#include <CGAL/Map_overlay_naive.h>
#include <CGAL/Boolean_operations_2.h>

#include <CGAL/leda_rational.h>
#include <LEDA/rat_window.h>
#include <CGAL/IO/Pm_Window_stream.h>
#include <CGAL/Bops_utility.h>

#if defined(LEDA_NAMESPACE)
using namespace leda;
#endif

typedef leda_rational                           NT;
typedef CGAL::Arr_leda_segment_exact_traits     Traits;

typedef Traits::Point_2                         Point;
typedef Traits::Curve_2                         Curve;
typedef Traits::X_curve                         X_curve;

typedef CGAL::Bop_default_dcel<Traits>          Dcel;
typedef CGAL::Planar_map_2<Dcel, Traits>        PM;

typedef CGAL::Map_overlay_default_notifier<PM>             
                                                Ovl_change_notification;
typedef CGAL::Map_overlay_2<PM, Ovl_change_notification>   MapOverlay;
typedef CGAL::Boolean_operations_2<MapOverlay>             Bops;

typedef CGAL::Pm_walk_along_line_point_location<PM>        PmWalkPL;

// global variables are used so that the redraw function for the LEDA window
// can be defined to draw information found in these variables.
static PmWalkPL pm_walk1, pm_walk2;
static PM pm1(&pm_walk1); 
static PM pm2(&pm_walk2);
static CGAL::Window_stream
  W(500, 500, "CGAL - Segment Boolean-Operations Demo");

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

/*void  draw_and_locate_maps (Bops& bops , 
                            const PM& pm1, 
                            const PM& pm2, 
                            CGAL::Window_stream& W)
{  
  W<<CGAL::BLUE;
  W<<pm1;
  W<<CGAL::RED;
  W<<pm2;

  for (;;) {
    double  x,y;

    int b = W.read_mouse(x,y);
    if (b==7) 
      break;
    
    W.set_node_width(3);
    W.clear();
    W<<CGAL::BLUE;
    W<<pm1;
    W<<CGAL::RED;
    W<<pm2;

    leda_color fg_col = leda_violet;
    W.set_color(fg_col);
    W.set_fill_color(fg_col);

    if (b == 4 || b == 5 || b == 6){
      Faces_container      faces_result;
      Halfedges_container  halfedges_result;
      Vertices_container   vertices_result;

      //Bops bops = Bops(pm1,pm2);
      
      if (b==4)
        bops.intersection(faces_result,halfedges_result,vertices_result);

      else if (b==5)       
        bops.Union(faces_result,halfedges_result,vertices_result);
      
      else if (b==6)
        bops.symmetric_difference(faces_result, halfedges_result,
                                  vertices_result);

      std::cout<<"number of faces is "<< faces_result.size()<<std::endl;
      std::cout<<"The number of resulting halfedges is "
               <<halfedges_result.size()<<std::endl;
      std::cout<<"The number of resulting vertices is "
               <<vertices_result.size()<<std::endl;

      for (list<PM::Face_const_handle>::iterator  f_iter =
           faces_result.begin(); 
           f_iter != faces_result.end(); ++f_iter){
        PM::Face_const_handle fh = *f_iter;
        if (fh->is_unbounded())
          continue;
        
        W.set_color(fg_col);
        W.set_fill_color(fg_col);
        
        PM::Ccb_halfedge_const_circulator cc = fh->outer_ccb();
        leda_list<leda_point> points_list;
        do {
          leda_point p = cc->source()->point().to_point();
          
          points_list.push_back(p);
        } while (++cc != fh->outer_ccb());
        W.draw_filled_polygon(points_list);
        
        for (PM::Holes_const_iterator hit = fh->holes_begin(); hit !=
             fh->holes_end(); ++hit) {
          leda_list<leda_point> holes_points_list;
          PM::Ccb_halfedge_const_circulator cc(*hit);
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
        } while (++cc != fh->outer_ccb());
        }
      
      for (list<PM::Halfedge_const_handle>::iterator h_iter =
           halfedges_result.begin(); 
           h_iter != halfedges_result.end(); ++h_iter, ++h_iter){
        PM::Halfedge_const_handle h = *h_iter;

        W.set_color(fg_col);
        W.set_node_width(4);
        W << h->curve();
      }
      
      for (list<PM::Vertex_const_handle>::iterator  v_iter =
           vertices_result.begin(); 
           v_iter != vertices_result.end(); v_iter++){
        W.set_node_width(4);
        
        PM::Vertex_const_handle vh = *v_iter;
        //Vertex_const_iterator vh = v_iter;
        
        W.set_color(leda_violet);
        W << vh->point();
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

  file.close();
  
  Traits traits;
  CGAL::sweep_to_construct_planar_map_2(curves.begin(), 
                                        curves.end(),
                                        traits, pm);
  W << pm;  
}

template <class Arrangement>
void calc_window_size(const Arrangement &arr, double &min_x, double &max_x,
                      double &min_y)
{
  Vertex_const_iterator v_iter = arr.vertices_begin();
  
  max_x = min_x = CGAL::to_double(v_iter->point().xcoordD());
  min_y = CGAL::to_double(v_iter->point().ycoordD());

  for (v_iter = arr.vertices_begin(); v_iter != arr.vertices_end(); v_iter++){
    NT x = v_iter->point().xcoordD(), y = v_iter->point().ycoordD();
    
    double  dx, dy;
    dx = CGAL::to_double(x);
    dy = CGAL::to_double(y);
    if (dx > max_x) max_x = dx;
    if (dx < min_x) min_x = dx;
    if (dy < min_y) min_y = dy;
  }
}*/

int main(int argc, char* argv[])
{
  double x0=-700,x1=700,y0=-700;
  
  if (argc != 3) {
    std::cout << "usage: Segment_sweep_ovl_from_file filename\n";
    exit(1);
  }
  
  CGAL::Bops_utility<PM,NT> utility;
  
  utility.scan_segment_planar_map(argv[1],pm1);
  utility.scan_segment_planar_map(argv[2],pm2);

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
      std::cout << "No edges were inserted to the first planar map. "
                << "First Planar map is empty. Exiting.";
      std::cout << std::endl;
    }
  
  // if second map is empty
  if (pm2.halfedges_begin() == pm2.halfedges_end()) 
    {
      std::cout << std::endl;
      std::cout << "No edges were inserted to the first planar map. "
                << "First Planar map is empty. Exiting.";
      std::cout << std::endl;
    }
  
  if (pm1.halfedges_begin() != pm1.halfedges_end() && 
      pm2.halfedges_begin() != pm2.halfedges_end() )
    utility.draw_and_locate_maps(bops,pm1,pm2,W);
 
  return 0;  
}

#endif  // CGAL_USE_LEDA

