// Constructs an arrangement of segments and circular arcs from a file.

// File format is:
// #number_of_arcs
// C x_centre y_centre squared_radius
//    or:
// A x_centre y_centre squared_radius x_source y_source x_target y_target
//    or:
// S x_source y_source x_target y_target

#include <CGAL/basic.h>

#ifndef CGAL_USE_LEDA
#include <iostream>
int main()
{
  std::cout << "Sorry, this demo needs LEDA for visualisation.";
  std::cout << std::endl;
  return 0;
}

#else

#include <CGAL/Cartesian.h>
#include <CGAL/Timer.h>
#include <CGAL/Arr_2_bases.h>
#include <CGAL/Arr_2_default_dcel.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Arr_conic_traits_2.h>
#include <CGAL/IO/Conic_arc_2_Window_stream.h>
#include <CGAL/IO/Window_stream.h>
#include <CORE/BigInt.h>
#include <CGAL/CORE_Expr.h>
#include <CGAL/Draw_preferences.h>

#include <fstream>

typedef CORE::BigInt                            CfNT;
typedef CGAL::Cartesian<CfNT>                   Int_kernel;
typedef Int_kernel::Point_2                     Int_point_2;
typedef Int_kernel::Circle_2                    Int_circle_2;
typedef Int_kernel::Segment_2                   Int_segment_2;

typedef CORE::Expr                              CoNT;
typedef CGAL::Cartesian<CoNT>                   Alg_kernel;

typedef CGAL::Arr_conic_traits_2<Int_kernel,
                                 Alg_kernel>    Traits_2;

typedef Traits_2::Point_2                       Point_2;
typedef Traits_2::Curve_2                       Curve_2;
typedef Traits_2::X_monotone_curve_2            X_monotone_curve_2;
typedef std::list<Curve_2>                      CurveList;

typedef CGAL::Pm_default_dcel<Traits_2>             Dcel;
typedef CGAL::Planar_map_2<Dcel,Traits_2>           Pm_2;
typedef CGAL::Planar_map_with_intersections_2<Pm_2> Pmwx_2;

// global variables are used so that the redraw function for the LEDA window
// can be defined to draw information found in these variables.
static Pmwx_2 arr;
static CGAL::Window_stream
W(400, 400, "CGAL - Segments and Circular Arcs Arrangement Demo");
 
CGAL_BEGIN_NAMESPACE
Window_stream& operator<<(Window_stream& os, Pmwx_2 &A)
{
  My_Arr_drawer<Pmwx_2, Pmwx_2::Ccb_halfedge_circulator, 
    Pmwx_2::Holes_iterator> drawer(os);
  
  draw_pm(arr, drawer, os);
  return os;
}
CGAL_END_NAMESPACE

// Redraw function for the LEDA window: used automatically when 
// the window reappears
void redraw(CGAL::Window_stream * wp) 
{ wp->start_buffering();
  wp->clear();
  // draw arragnement
  *wp << arr;
  wp->flush_buffer();
  wp->stop_buffering();
}

// The main:
int main(int argc, char* argv[])
{
  CGAL::Timer insrt_t;

  if (argc<2) {
    std::cerr << "usage: Circle_arr_from_file <filename>" << std::endl;
    exit(1);
  }

  // Read the number of arcs from a file.
  int  n_arcs;
  std::ifstream f(argv[1]);
  f >> n_arcs;

  double min_x = -1, max_x=1;   // For adjusting the window size.
  double min_y = -1, max_y = 1; // For adjusting the window size.
  char   type;
  int    i_arc;
  
  for (i_arc = 0; i_arc < n_arcs; i_arc++)
  {
    // Read the arc type.
    f >> type;

    std::cout << "Inserting arc no. " << i_arc + 1;

    // A full circle (c) or a circular arc (a):
    if (type == 'c' || type == 'C' || type == 'a' || type == 'A')
    {
      // Read the circle, using the format "x0 y0 r^2"
      int    x0, y0, r2;
    
      f >> x0 >> y0 >> r2;

      Int_point_2    center = Int_point_2 (CfNT(x0), CfNT(y0));
      Int_circle_2   circle = Int_circle_2 (center, CfNT(r2));

      if (type == 'c' || type == 'C')
      {
	std::cout << " (full circle)." << std::endl;

	insrt_t.start();
	arr.insert (Curve_2(circle));
	insrt_t.stop();
      }
      else
      {
	std::cout << " (circular arc)." << std::endl;

	// Read the end points.
	int    x1, y1, x2, y2;

	f >> x1 >> y1 >> x2 >> y2;

	Point_2      source = Point_2 (CoNT(x1), CoNT(y1));
	Point_2      target = Point_2 (CoNT(x2), CoNT(y2));

	insrt_t.start();
	arr.insert (Curve_2 (circle, CGAL::CLOCKWISE, source, target));
	insrt_t.stop();
      }

      // Check whether we need to resize the screen.
      double dx = x0;
      double dy = y0;
      double dr = r2;

      if (min_x > dx - dr) 
	min_x = dx - dr;
      if (max_x < dx + dr)
	max_x = dx + dr;
      if (min_y > dy - dr)
	min_y = dy - dr;
      if (max_y < dy + dr)
	max_y = dy + dr;
    }
    // A segment (s):
    else if (type == 's' || type == 'S')
    {
      std::cout << " (segment)." << std::endl;
      
      // Read the end points.
      int    x1, y1, x2, y2;

      f >> x1 >> y1 >> x2 >> y2;

      Int_point_2      source = Int_point_2 (CfNT(x1), CfNT(y1));
      Int_point_2      target = Int_point_2 (CfNT(x2), CfNT(y2));

      insrt_t.start();
      arr.insert (Curve_2 (Int_segment_2 (source, target)));
      insrt_t.stop();

      // Check whether we need to resize the screen.
      double xmin = CGAL::to_double(x1 < x2 ? x1 : x2);
      double xmax = CGAL::to_double(x1 > x2 ? x1 : x2);
      double ymin = CGAL::to_double(y1 < y2 ? y1 : y2);
      double ymax = CGAL::to_double(y1 > y2 ? y1 : y2);

      if (min_x > xmin) 
	min_x = xmin;
      if (max_x < xmax)
	max_x = xmax;
      if (min_y > ymin)
	min_y = ymin;
      if (max_y < ymax)
	max_y = ymax;
    }
    else
    {
      std::cout << "Unknown arc type '" << type 
		<< "' - Stopping here." << std::endl;
      exit(1);
    }
  }
  f.close();

  // Draw the arrangement.
  double  diff_x = max_x - min_x;
  double  diff_y = max_y - min_y;
  double  max_diff = (diff_x > diff_y) ? diff_x : diff_y;
  double  margin_percentage = 8;
  double  margin = (max_diff * margin_percentage) / 100.0;

  W.init (min_x - margin, max_x + max_diff + 2*margin,
	  min_y - margin);

  W.set_redraw(redraw);
  W.set_mode(CGAL_LEDA_SCOPE::src_mode);
  W.set_node_width(3);
  W.button("finish",10);
  W.open_status_window();
  W.display();

  W << arr;
  std::cout << "Total insertion time: " <<  insrt_t.time() << std::endl;

  //POINT LOCATION
  W.set_status_string( "Left mouse button - query point." );
  Point_2 p;

  Pmwx_2::Halfedge_handle e;
  
  for (;;) {
    double x,y;
    int b=W.read_mouse(x,y);
    if (b==10) break;
    else
      p=Point_2(x,y);

    W << arr;
    
    Pmwx_2::Locate_type lt;
    e = arr.locate(p,lt);

    Pmwx_2::Face_handle fh=e->face();
    //Pmwx_2::Ccb_halfedge_circulator cc(e);
    Pmwx_2::Ccb_halfedge_circulator cc;

    if (fh != arr.unbounded_face()) {
      cc=fh->halfedge_on_outer_ccb();
      do {
        W << cc->curve();
      } while (++cc != fh->halfedge_on_outer_ccb());
    }

    Pmwx_2::Holes_iterator hit=fh->holes_begin(), eit=fh->holes_end();
    for (;hit!=eit;++hit) {
      cc=*hit;
      do {
        W << cc->curve();
      } while (++cc != *hit);
    }
      
  }

  return 0;  
}

#endif
