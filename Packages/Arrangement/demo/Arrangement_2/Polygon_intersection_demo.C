#include <iostream>

#ifndef CGAL_USE_LEDA
int main(int argc, char* argv[])
{

  std::cout << "Sorry, this demo needs LEDA for visualisation.";
  std::cout << std::endl;

  return 0;
}
 
#else

// Uses types associated with Polygon as declared in the interface of
// Polygon_intersection.
#include "Polygon_intersection.h"

#include <CGAL/IO/Window_stream.h>
#include <CGAL/IO/leda_window.h>
#include <LEDA/polygon.h>
#include <LEDA/list.h>
#include <LEDA/string.h>

typedef CGAL::Window_stream Window;

// A convertion function so we can output polygons to the screen.
leda_polygon cgal_to_leda_polygon(Polygon& poly)
{
  leda_list<leda_point> plist;
  leda_polygon leda_poly;

  for (Polygon::Vertex_iterator it=poly.vertices_begin();
       it != poly.vertices_end();
       it++) {
    plist.push_back(leda_point(CGAL::to_double(it->x()), 
			       CGAL::to_double(it->y())));
  }
  
  return leda_polygon(plist);
}

// read a polygon from mouse 
void PolygonRead(Polygon &poly, Window &W)
{
  poly.erase(poly.vertices_begin(), poly.vertices_end());
  leda_list<leda_point> leda_poly = W.read_polygon();
  leda_list<leda_point>::iterator it;
  for (it = leda_poly.begin(); it != leda_poly.end(); it++)
    {
      poly.push_back(Point((*it).xcoord(), (*it).ycoord()));
    }
}  
 
// A redraw function
void redraw(Window &W, Polygon &poly1, Polygon &poly2,
            Polygon_list &poly_list)
{
  W.clear();

  W << CGAL::PURPLE;
  for (Polygon_list::iterator pit = poly_list.begin(); 
       pit != poly_list.end(); 
       pit++) {
    W.draw_filled_polygon(cgal_to_leda_polygon(*pit));
  }
  W.set_line_width(2);
  W << CGAL::RED << poly1;
  W << CGAL::BLUE << poly2;
}

int main(int argc, char* argv[])
{
  Polygon poly1, poly2;
  Polygon_list in_poly_list;
  Polygon_list out_poly_list;
  Window W(600, 600);

  // Define constants for buttons
  enum { CLEAR, POLY1, POLY2, INTERSECT, REFRESH, QUIT};

  // Set window parameters
  W.init(0, 1000, 0);
  W.set_mode(leda_src_mode);
  W.set_node_width(3);
    
  // Define Buttons
  W.button("Clear",        CLEAR);
  W.button("Insert Poly1", POLY1);
  W.button("Insert Poly2", POLY2);
  W.button("Intersect ",   INTERSECT);
  W.button("Refresh",      REFRESH);
  W.button("Quit",         QUIT);
  W.display();

  leda_panel welcome_panel("Welcome");

  welcome_panel.text_item(
	     "This demo calculates the intersection of the interiors of two "
	     "simple polygons. "
	     "Insert each polygon by pressing [Insert Poly1] "
	     "or [Insert Poly2]. "
	     "Insert points by clicking the left mouse button. "
	     "Close the polygon by clicking the right mouse button.");
  welcome_panel.button("continue");

  welcome_panel.open(W);
  W.clear();

  // Main demo loop
  double x,y; 
  while( true )
    {
      int b = W.read_mouse(x,y);
      // Quit button
      if (b == QUIT) break;

      // Clear button
      if (b == CLEAR)
        {
	  poly1.erase(poly1.vertices_begin(), poly1.vertices_end());
	  poly2.erase(poly2.vertices_begin(), poly2.vertices_end());
	  out_poly_list.clear();
        }

      // Inserting polygon 1
      if (b == POLY1)
        {
	  out_poly_list.clear();
	  poly1.erase(poly1.vertices_begin(), poly1.vertices_end());
	  redraw(W, poly1, poly2, out_poly_list);

	  W << CGAL::RED;
	  PolygonRead(poly1, W);
	  out_poly_list.clear();
        }

      // Inserting polygon 2
      if (b == POLY2)
        {
	  out_poly_list.clear();
	  poly2.erase(poly2.vertices_begin(), poly2.vertices_end());
	  redraw(W, poly1, poly2, out_poly_list);

	  W << CGAL::BLUE;
	  PolygonRead(poly2, W);
	  out_poly_list.clear();
        }

      // Polygon intersection
      if (b == INTERSECT)
        {
	  out_poly_list.clear();
	  in_poly_list.clear();
	  in_poly_list.push_back(poly1);
	  in_poly_list.push_back(poly2);
	  intersect_polygons(in_poly_list, out_poly_list);
        }

      // Perform redraw for whatever operation was performed.
      // In particular, for a refresh operation.
      redraw(W, poly1, poly2, out_poly_list);
    } 
  return 0;
}

#endif
