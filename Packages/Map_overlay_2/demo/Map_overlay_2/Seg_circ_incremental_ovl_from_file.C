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

#include <CGAL/Pm_walk_along_line_point_location.h>
#include <CGAL/Pm_with_intersections.h>

#include <CGAL/Arr_segment_circle_traits.h>
#include <CGAL/IO/Segment_circle_Window_stream.h>

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

typedef leda_real                             NT;
typedef CGAL::Arr_segment_circle_traits<NT>   Traits; 

typedef Traits::Point                                  Point;
typedef Traits::Segment                                Segment;
typedef Traits::Circle                                 Circle;
typedef Traits::Conic                                  Conic;
typedef Traits::Curve                                  Curve; 
typedef Traits::X_curve                                X_curve;

typedef CGAL::Map_overlay_default_dcel<Traits>                Dcel;
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

/*void  draw_and_locate_maps (const Pmwx& pmwx, 
                            const Pmwx& pmwx1, 
                            const Pmwx& pmwx2, 
                            CGAL::Window_stream& W)
{ 
  for (;;) {
    double  x,y;

    int b = W.read_mouse(x,y);
    if (b == 4) 
      break;

    W.clear();
    W.set_line_width(1);
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

      //if (f1 == NULL || f2 == NULL)
      //  std::cout<<"NULL faces pointers\n";
 
      leda_color fg_col;
      if (f1 != NULL && !(f1->is_unbounded()) && f2 != NULL && !(f2->is_unbounded()) )
        fg_col = leda_violet;
      else if (f1 != NULL && !(f1->is_unbounded()))
        fg_col = leda_blue;
      else if (f2 != NULL && !(f2->is_unbounded()))
        fg_col = leda_red;
      else
        fg_col = leda_orange;
      
      W.set_color(fg_col);
      //W.set_fill_color(fg_col);

      Pmwx::Ccb_halfedge_const_circulator cc = fh->outer_ccb();
      W.set_line_width(3);
      do {
        W << (*cc).curve();
        W << cc->source()->point();
      } while (++cc != fh->outer_ccb());
      
      for (Pmwx::Holes_const_iterator hit = fh->holes_begin(); hit != fh->holes_end(); ++hit, ++hit) {
        Pmwx::Ccb_halfedge_const_circulator cc(*hit);
        do{
          W << cc->curve();
        } while (++cc != *hit);
      }
      
      for (Pmwx::Vertex_const_iterator  v_iter = pmwx.vertices_begin(); 
           v_iter != pmwx.vertices_end(); v_iter++){
        if (v_iter->get_first_vertex_above() != NULL && v_iter->get_second_vertex_above() != NULL)
          W.set_color(leda_violet);
        else if (v_iter->get_first_vertex_above() != NULL && v_iter->get_second_vertex_above() == NULL)
          W << CGAL::BLUE;
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

void  scan_pmwx(const char* filename, Pmwx& pmwx)
{
  std::ifstream file(filename);
  int num_curves;
  file >> num_curves;
  char   type;
  int    i_arc;
  
  for (i_arc = 0; i_arc < num_curves; ++i_arc)
  {
    // Read the arc type.
    file >> type;

    std::cout << "Inserting arc no. " << i_arc + 1;

    // A full circle (c) or a circular arc (a):
    if (type == 'c' || type == 'C' || type == 'a' || type == 'A')
      {
      // Read the circle, using the format "x0 y0 r^2"
      NT  x0, y0, r2;
    
      file >> x0 >> y0 >> r2;
    
      Circle      circle (Point (x0, y0), r2, CGAL::CLOCKWISE);

      if (type == 'c' || type == 'C')
      {
	std::cout << " (full circle)." << std::endl;
	pmwx.insert (Curve(circle));
      }
      else
      {
	std::cout << " (circular arc)." << std::endl;

	// Read the end points.
	NT  x1, y1, x2, y2;

	file >> x1 >> y1 >> x2 >> y2;

	Point      source (x1, y1);
	Point      target (x2, y2);
        
	pmwx.insert (Curve (circle, source, target));
      }
    }
    // A segment (s):
    else if (type == 's' || type == 'S')
    {
      std::cout << " (segment)." << std::endl;
      // Read the end points.
      NT  x1, y1, x2, y2;

      file >> x1 >> y1 >> x2 >> y2;

      Point      source (x1, y1);
      Point      target (x2, y2);
      
      pmwx.insert (Curve (Segment (source, target)));
    }
    else
      {
      std::cout << "Unknown arc type '" << type 
		<< "' - Stopping here." << std::endl;
      exit(1);
    }
  }
  
  file.close();
  
  W << pmwx;  
}*/

int main(int argc, char* argv[])
{
  double x0=-200,x1=200,y0=-200;
  
  if (argc != 3) {
    std::cout << "usage: Segment_sweep_ovl_from_file filename\n";
    exit(1);
  }

  CGAL::Ovl_utility<Pmwx,NT> utility;
  
  utility.scan_seg_circ_pmwx(argv[1],pmwx1);
  utility.scan_seg_circ_pmwx(argv[2],pmwx2);
 
  utility.calc_window_size(pmwx1,x0,x1,y0);
  utility.calc_window_size(pmwx2,x0,x1,y0);
    
  W.init(x0,x1,y0);
  W.set_redraw(redraw);
  W.set_mode(leda_src_mode);
  W.set_node_width(3);
  W.open_status_window();
  W.button("Exit",4);
  W.display();
  
  W << CGAL::RED;
  W << pmwx1;
  W << CGAL::BLUE;
  W << pmwx2;
  
  MapOverlay_incremental   ovl_incremental, pmwx_incremental1, pmwx_incremental2;
  MapOverlay map1(pmwx1, &pmwx_incremental1);
  MapOverlay map2(pmwx2, &pmwx_incremental2);
  MapOverlay map_overlay(map1, map2, &ovl_incremental);
  
  std::cout<<"Locate Overlay Face:"<<endl;
  std::cout<<"Purple Face - an overlay face laying under two bounded faces"<<std::endl;
  std::cout<<"Blue Face - an overlay face laying under a bounded face of the first map and the unbounded face of the second"<<std::endl;
  std::cout<<"Red Face - an overlay face laying under the unbounded face of the first map and a bounded face of the second"<<std::endl;
  std::cout<<"Orange Face - an overlay face laying under the unbounded faces of both maps"<<endl;

  // Point Location Queries
  W.set_status_string("Map Overlay. Enter a point with left button.");
  
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
