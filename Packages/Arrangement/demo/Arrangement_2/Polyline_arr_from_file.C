// demo/Arrangement_2/Polyline_arr_from_file.C
//
// constructs an arrangement of polylines from file
// We use the leda traits (therefore we use leda functions).

#include "short_names.h"
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
#include <CGAL/MP_Float.h>
#include <CGAL/Quotient.h>
#include <CGAL/Arr_2_bases.h>
#include <CGAL/Arr_2_default_dcel.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Arr_segment_cached_traits_2.h>
#include <CGAL/Arr_polyline_traits_2.h>

typedef CGAL::Quotient<CGAL::MP_Float>                NT;
typedef CGAL::Cartesian<NT>                           Kernel;
typedef Kernel::Segment_2                             Segment_2;
typedef CGAL::Arr_segment_cached_traits_2<Kernel>     Seg_traits;
typedef CGAL::Arr_polyline_traits_2<Seg_traits>       Traits;

typedef Traits::Point_2                               Point_2;
typedef Traits::Curve_2                               Curve_2;
typedef Traits::X_monotone_curve_2                    X_monotone_curve_2;

typedef CGAL::Pm_default_dcel<Traits>                 Dcel;
typedef CGAL::Planar_map_2<Dcel,Traits>               Pm_2;
typedef CGAL::Planar_map_with_intersections_2<Pm_2>   Pmwx_2;

Pmwx_2 pmwx;           // The arrangement, defined as a global variable.

#include <CGAL/IO/Window_stream.h>

// draw a polyline, with points as 'x's
CGAL::Window_stream& operator<< (CGAL::Window_stream& os, 
				 const X_monotone_curve_2& cv)
{
  X_monotone_curve_2::const_iterator ps = cv.begin();
  X_monotone_curve_2::const_iterator pt = ps; pt++;

  while (pt != cv.end())
  {
    os << Segment_2(*ps, *pt);
    ps++; pt++;
  }

  return (os);
}

CGAL::Window_stream& operator<< (CGAL::Window_stream& os, Pmwx_2& _pmwx)
{
  Pmwx_2::Edge_iterator ei;
  os << CGAL::BLUE;
  for (ei = _pmwx.edges_begin(); ei != _pmwx.edges_end(); ++ei)
    os << (*ei).curve();

  Pmwx_2::Vertex_iterator vi;
  os << CGAL::RED;
  for (vi = _pmwx.vertices_begin(); vi != _pmwx.vertices_end(); ++vi)
    os << (*vi).point();

  return (os);
}

// Read the curves from a file and compute their bounding box.
void read_curves (const char* filename,
		  std::list<Curve_2>& curves,
		  CGAL::Bbox_2& bbox)
{  
  std::ifstream file(filename);
  curves.clear();

  int                num_polylines, num_segments;
  NT                 x, y;
  std::list<Point_2> points;
  int                i, j;

  file >> num_polylines;
  for (i = 0; i < num_polylines; i++) 
  {
    file >> num_segments;
    points.clear();
    for (j = 0; j < num_segments; j++)
    {
      file >> x >> y;
      points.push_back (Point_2(x,y));
    }

    Curve_2   polyline(points.begin(), points.end());
    curves.push_back(polyline);

    if (i == 0)
      bbox = polyline.bbox();
    else
      bbox = bbox + polyline.bbox();
  }

  return;
}

// redraw function for LEDA window. used automatically when window reappears
void redraw(CGAL::Window_stream * wp) 
{
  wp->start_buffering();
  wp->clear();
  // draw arragnement
  *wp << pmwx;
  wp->flush_buffer();
  wp->stop_buffering();
}

int main(int argc, char* argv[])
{
  if (argc != 2) {
    std::cout << "usage: Polyline_arr_from_file <filename>" << std::endl;
    exit(1);
  }

  // Read the polyline curves from the input file.
  std::list<Curve_2>                 curves;
  std::list<Curve_2>::const_iterator cv_it;
  CGAL::Bbox_2                       bbox;

  read_curves (argv[1],
	       curves,
	       bbox);

  for (cv_it = curves.begin(); cv_it != curves.end(); cv_it++)
    pmwx.insert (*cv_it);

  // Initialize the display window.
  float x_range = bbox.xmax() - bbox.xmin();
  float y_range = bbox.ymax() - bbox.ymin();
  float width = 640;
  float height = (y_range * width) / x_range;
    
  CGAL::Window_stream W (static_cast<int>(width),
			 static_cast<int>(height),
			 "CGAL - Polyline Arrangement Demo");

  float min_range = (x_range < y_range) ? x_range : y_range;
  float x_margin = min_range / 5;
  float y_margin = (height * x_margin) / width;
        
  float x0 = bbox.xmin() - x_margin;
  float x1 = bbox.xmax() + x_margin;
  float y0 = bbox.ymin() - y_margin;
  W.init(x0, x1, y0);   // logical window size 

  W.set_redraw(redraw);
  W.set_mode(CGAL_LEDA_SCOPE::src_mode);
  W.set_node_width(3);
  const int THE_BUTTON = 10;
  W.button("  Quit  ", THE_BUTTON);
  W.open_status_window();
  W.display(leda_window::center,leda_window::center);

  // Draw the arrangement.
  W << pmwx;
  
  // Point Location part.
  Pmwx_2::Halfedge_handle e;
  double  x,y;
  Point_2 pnt;

  W.set_status_string("Enter a point with left button. Press Quit to quit.");  
  W.set_button_label(THE_BUTTON, "  Quit  ");
  while (W.read_mouse(x,y) != THE_BUTTON) 
  {
    // Read the query point from the mouse.
    pnt = Point_2(x,y);
    W << pmwx;
    
    Pmwx_2::Locate_type lt;
    e = pmwx.locate(pnt ,lt);

    // Color the face containing the query point on the screen.
    W << CGAL::GREEN;
    Pmwx_2::Face_handle f=e->face();
    if (f->does_outer_ccb_exist()) 
    {
      Pmwx_2::Ccb_halfedge_circulator cc = f->outer_ccb();
      do 
      {
	W << cc->curve();
      } while (++cc != f->outer_ccb());
    }

    for (Pmwx_2::Holes_iterator ho = f->holes_begin(), hoe = f->holes_end();
	 ho != hoe; ++ho) 
    {
      Pmwx_2::Ccb_halfedge_circulator cc = *ho; 
      do {
	W << cc->curve();
      } while (++cc != *ho);
    }
  }

  return (0);  
}

#endif
