// demo/Arrangement_2/Polyline_arr_from_mouse.C
//
//constructs an arrangement of polylines from user input

#include "short_names.h"

//constructs a ployline arrangement from CGAL window.
// We use the leda traits (therefore we use leda functions).
#include <CGAL/basic.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Arr_2_bases.h>
#include <CGAL/Arr_2_default_dcel.h>
#include <CGAL/Arr_polyline_traits.h>

#ifndef CGAL_USE_LEDA
int main()
{

  std::cout << "Sorry, this demo needs LEDA for visualisation.";
  std::cout << std::endl;

  return 0;
}

#else

#include <CGAL/leda_real.h>
#include <LEDA/string.h>

#include <CGAL/Draw_preferences.h>

typedef leda_real                                       NT;
typedef CGAL::Cartesian<NT>                             Kernel;
typedef CGAL::Arr_polyline_traits<Kernel>               Traits;

typedef Traits::Point_2                                 Point;
typedef Traits::X_curve_2                               X_curve;
typedef Traits::Curve_2                                 Curve;

typedef CGAL::Arr_base_node<X_curve>                    Base_node;
typedef CGAL::Arr_2_default_dcel<Traits>                Dcel;

typedef CGAL::Arrangement_2<Dcel,Traits,Base_node >     Arr_2;

// global variables are used so that the redraw function for the LEDA window
// can be defined to draw information found in these variables.
static Arr_2 Arr; 
static leda_string         text("");
static Curve               cv1;

leda_point to_leda_pnt(Point p)
{
  return leda_point(p.x().to_double(), p.y().to_double());
}

CGAL_BEGIN_NAMESPACE
// drawing functions
Window_stream & operator<<(Window_stream & os, const  Point & p)
{
  // conversion to leda_point in order to show it on screen
  return os << to_leda_pnt(p);
}

// draw a polyline, with points as 'x's
Window_stream & operator<<(Window_stream & os, const X_curve & c)
{
  X_curve::const_iterator sit, tit;
  sit = c.begin();
  tit = sit; tit++;
  os.draw_point(to_leda_pnt(*sit), leda_green); // draw first point
  for (; tit != c.end(); tit++, sit++) {
    // conversion to screen drawble segment
    os << leda_segment(to_leda_pnt(*sit), to_leda_pnt(*tit));
    os.draw_point(to_leda_pnt(*tit), leda_green);
  }
    
  return os;
}

Window_stream & operator<<(Window_stream & os, Arr_2 & arr)
{
  My_Arr_drawer< Arr_2,
                 Arr_2::Ccb_halfedge_circulator, 
                 Arr_2::Holes_iterator> drawer(os);
  draw_pm(arr, drawer, os);
  return os;
}
CGAL_END_NAMESPACE

void show_welcome_message(CGAL::Window_stream & os)
{
  text += "\\black ";
  text += "Click left button for polyline points.\\n ";
  text += "Click right button as last point in polygon.\\n ";
  text += " \\n ";
  text +=
    "Clicking close to a point, assumes the location is at the point.\\n ";
  text += "Lonely points will be discarded.\\n ";
 
  os.set_status_string("Press Begin to enter polylines.");
  os.redraw();
  // wait for button
}

// redraw function for LEDA window. used automatically when window reappears
void redraw(CGAL::Window_stream * wp) 
{
  wp->start_buffering();
  wp->clear();
  // display message if one exists
  if (text) wp->text_box(text);
  // draw currently inserted polyline
    *wp << CGAL::BLACK;
  if (cv1.size() > 1) *wp << cv1;
  if (cv1.size() == 1) *wp << cv1[0];
  // draw arragnement
  *wp << Arr;
  wp->flush_buffer();
  wp->stop_buffering();
}

int main()
{
  // window settings
  double x0=-400,x1=400,y0=-400;
  enum { THE_BUTTON = 10 };
  CGAL::Window_stream W(400, 400, "CGAL - 2D Polyline Arrangement Demo");
  W.init(x0,x1,y0);
  W.set_redraw(redraw);
  W.set_mode(leda_src_mode);
  W.set_node_width(3);
  W.button("  Begin  ",THE_BUTTON);
  W.open_status_window();
  W.display(leda_window::center,leda_window::center);

  // (1) welcome part
  show_welcome_message(W);
  while (W.read_mouse() != THE_BUTTON);
  text = ""; // empty message for redraw
  W.clear();
  W.set_status_string("Left - points. "
                      "Right - last point. Continue - point location");

  // (2) polylines insertion part
  Point  pnt, last_pnt;
  bool   should_exit = false,
         first_point = true;
  double x, y; // doubles are used to read off screen, but not in computations
  int    b; // button

  W.set_button_label(THE_BUTTON, "Continue");

  while (! should_exit) {
    b = W.read_mouse(x,y);
    last_pnt = pnt;
    pnt = Point(x,y);
  
    if (b == THE_BUTTON)
      should_exit = true;
    else if (b == MOUSE_BUTTON(1) || MOUSE_BUTTON(3)) {
      // looking for points in the vicinity of the click
      bool vicinity_point = false,
        last_in_polyline = false;      

      // first, looking in the points of this polyline
      Arr_2::Curve::iterator cit = cv1.begin();
      for(; ! vicinity_point && cit != cv1.end(); cit++) {
	if ( CGAL::squared_distance(pnt, *cit) < ((x1-x0)/50)*((x1-x0)/50) ) {
          pnt =* cit;
          vicinity_point = true;
        }
      }
      // check if last point was re-clicked
      if ( vicinity_point && last_pnt == pnt ) {
        last_in_polyline = true;
      }
      // second, looking in other points of the arrangement
      for(Arr_2::Vertex_iterator vi = Arr.vertices_begin();
	  ! vicinity_point && vi != Arr.vertices_end(); ++vi) {
	if ( CGAL::squared_distance(pnt, vi->point()) <
             ((x1-x0)/50)*((x1-x0)/50) )
        {
          pnt = vi->point(); 
          vicinity_point = true;
        }
      }
      // if last point was re-clicked ignore to avoid invalidity of polyline.
      if ( ! last_in_polyline )
	cv1.push_back(pnt);
      W << CGAL::BLACK;
      W << pnt;
      if ( ! first_point ) {
        W << leda_segment(last_pnt.x().to_double(), last_pnt.y().to_double(),
                          pnt.x().to_double(), pnt.y().to_double());
      }
      first_point = false;
    }
    // end of polyline 
    // on right click, or on "Continue"
    if (b == MOUSE_BUTTON(3) || (should_exit && cv1.size() > 1)) {
      if ( cv1.size() > 1) Arr.insert(cv1);
      //(at least 2 points, otherwise ignore)
      cv1.clear();
      W << Arr;
      first_point = true;
    }
  }

  W.redraw();

  // (3) Point Location part
  Arr_2::Halfedge_handle e;
  bool map_is_empty = false;

  if (Arr.halfedges_begin() == Arr.halfedges_end()) 
    map_is_empty = true;

  if (map_is_empty) {
    W.set_status_string("Arrangement is empty. Press Quit to quit.");
  }
  else {
    W.set_status_string("Enter a point with left button. Press Quit to quit.");
  }
  W.set_button_label(THE_BUTTON, "  Quit  ");
  //enable_button(THE_BUTTON);
  while (W.read_mouse(x,y) != THE_BUTTON) {
    pnt = Point(x, y);
    W << Arr;
    
    if ( ! map_is_empty ) {
      Arr_2::Locate_type lt;
      e = Arr.locate(pnt ,lt);

      // color the face on the screen
      W << CGAL::GREEN;
      Arr_2::Face_handle f=e->face();
      if (f->does_outer_ccb_exist()) {
	Arr_2::Ccb_halfedge_circulator cc=f->outer_ccb();
	do {
	  W << cc->curve();
	} while (++cc != f->outer_ccb());
      }
      for (Arr_2::Holes_iterator ho=f->holes_begin(),hoe = f->holes_end();
	   ho != hoe; ++ho) {
	Arr_2::Ccb_halfedge_circulator cc = *ho; 
	do {
	  W << cc->curve();
	} while (++cc != *ho);
      }
    } // if ! map_is_empty
  } // while

  return 0;  
}

#endif
