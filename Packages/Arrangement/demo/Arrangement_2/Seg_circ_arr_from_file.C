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
int main(int argc, char* argv[])
{

  std::cout << "Sorry, this demo needs LEDA for visualisation.";
  std::cout << std::endl;

  return 0;
}

#else

#include <CGAL/Cartesian.h>
#include <fstream>
#include <CGAL/Timer.h>

#include <CGAL/Arr_segment_circle_traits.h>
#include <CGAL/IO/Segment_circle_Window_stream.h>

#include <CGAL/Arr_2_bases.h>
#include <CGAL/Arr_2_default_dcel.h>
#include <CGAL/Arrangement_2.h>

#include <CGAL/leda_real.h>
#include <CGAL/Draw_preferences.h>

typedef CGAL::Arr_segment_circle_traits<leda_real> Traits; 

typedef Traits::Point                                  Point;
typedef Traits::Segment                                Segment;
typedef Traits::Circle                                 Circle;
typedef Traits::Conic                                  Conic;
typedef Traits::Curve                                  Curve; 
typedef Traits::X_curve                                X_curve;

typedef CGAL::Arr_base_node<Curve>                     Base_node;
typedef CGAL::Arr_2_default_dcel<Traits>               Dcel;
typedef CGAL::Arrangement_2<Dcel,Traits,Base_node >    Arr_2;

// global variables are used so that the redraw function for the LEDA window
// can be defined to draw information found in these variables.
static Arr_2               arr;
static CGAL::Window_stream W(400, 400, 
			    "CGAL - Segments and Circular Arcs Arrangement Demo");
 
CGAL_BEGIN_NAMESPACE
Window_stream& operator<<(Window_stream& os, Arr_2 &A)
{
  My_Arr_drawer< Arr_2,
                 Arr_2::Ccb_halfedge_circulator, 
                 Arr_2::Holes_iterator> drawer(os);
  
  draw_pm(arr, drawer, os);
  
  return os;
}
CGAL_END_NAMESPACE

// Redraw function for the LEDA window: used automatically when window reappears
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
      leda_real  x0, y0, r2;
    
      f >> x0 >> y0 >> r2;
    
      Circle      circle (Point (x0, y0), r2);

      if (type == 'c' || type == 'C')
      {
	std::cout << " (full circle)." << std::endl;

	insrt_t.start();
	arr.insert (Curve(circle));
	insrt_t.stop();

      }
      else
      {
	std::cout << " (circular arc)." << std::endl;

	// Read the end points.
	leda_real  x1, y1, x2, y2;

	f >> x1 >> y1 >> x2 >> y2;

	Point      source (x1, y1);
	Point      target (x2, y2);

	insrt_t.start();
	arr.insert (Curve (circle, source, target));
	insrt_t.stop();
      }

      // Check whether we need to resize the screen.
      double dx = CGAL::to_double(x0);
      double dy = CGAL::to_double(y0);
      double dr = sqrt(CGAL::to_double(r2));

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
      leda_real  x1, y1, x2, y2;

      f >> x1 >> y1 >> x2 >> y2;

      Point      source (x1, y1);
      Point      target (x2, y2);
 
      insrt_t.start();
      arr.insert (Curve (Segment (source, target)));
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
  if (max_y - min_y < max_x - max_y)
    W.init (min_x - 1, max_x + 1, min_y - 1);
  else
    W.init (min_x - 1, min_x + (max_y - min_y) +  1, min_y - 1);

  W.set_redraw(redraw);
  W.set_mode(leda_src_mode);
  W.set_node_width(3);
  W.button("finish",10);
  W.open_status_window();
  W.display();

  W << arr;
  std::cout << "Total insertion time: " <<  insrt_t.time() << std::endl;

  //POINT LOCATION
  W.set_status_string( "Left mouse button - query point." );
  Point p;

  Arr_2::Halfedge_handle e;
  
  for (;;) {
    double x,y;
    int b=W.read_mouse(x,y);
    if (b==10) break;
    else
      p=Point(x,y);

    W << arr;
    
    Arr_2::Locate_type lt;
    e = arr.locate(p,lt);

    Arr_2::Face_handle fh=e->face();
    //Arr_2::Ccb_halfedge_circulator cc(e);
    Arr_2::Ccb_halfedge_circulator cc;

    if (fh != arr.unbounded_face()) {
      cc=fh->halfedge_on_outer_ccb();
      do {
        W << cc->curve();
      } while (++cc != fh->halfedge_on_outer_ccb());
    }

    Arr_2::Holes_iterator hit=fh->holes_begin(), eit=fh->holes_end();
    for (;hit!=eit;++hit) {
      cc=*hit;
      do {
        W << cc->curve();
      } while (++cc != *hit);
    }
      
  }

  return 0;  
}

#endif // CGAL_USE_LEDA
