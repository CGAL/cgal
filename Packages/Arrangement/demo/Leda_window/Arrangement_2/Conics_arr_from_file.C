#include <CGAL/basic.h>

#ifndef CGAL_USE_LEDA
#include <iostream>
int main()
{
  std::cout << "Sorry, this demo needs LEDA." << std::endl;
  return 0;
}

#else

#include <CGAL/Cartesian.h>
#include <CGAL/leda_real.h>
#include <CGAL/Arr_conic_traits_2.h>
#include <CGAL/Timer.h>

#include <CGAL/Pm_default_dcel.h>
#include <CGAL/Planar_map_2.h>
#include <CGAL/Pm_with_intersections.h>
#include <CGAL/IO/Window_stream.h>
#include <CGAL/IO/Pm_iostream.h>
#include <CGAL/IO/Pm_Window_stream.h>
#include <CGAL/IO/Conic_arc_2_Window_stream.h>

#include <CGAL/Pm_default_point_location.h>
#include <CGAL/Pm_walk_along_line_point_location.h>
#include <CGAL/Pm_naive_point_location.h>
#include <CGAL/Pm_dummy_point_location.h>

#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <list>

#include <CGAL/Draw_preferences.h>

enum FormatId  {
  FORMAT_RAT = 0,
  FORMAT_INT,
  FORMAT_FLT,
  FORMAT_R,
  FORMAT_I,
  FORMAT_F
};

typedef leda_real                                       NT;
typedef CGAL::Cartesian<NT>                             Kernel;
typedef CGAL::Arr_conic_traits_2<Kernel>                Traits;
typedef CGAL::Pm_default_dcel<Traits>                   Dcel;
typedef CGAL::Planar_map_2<Dcel,Traits>                 Pm;
typedef CGAL::Planar_map_with_intersections_2<Pm>       Pmwx;

typedef Traits::Point_2                                 Point_2;
typedef Traits::Curve_2                                 Curve_2;
typedef Traits::Circle_2                                Circle_2;
typedef Traits::Segment_2                               Segment_2;
typedef Traits::X_monotone_curve_2                      X_monotone_curve_2;
typedef std::list<Curve_2>                              CurveList;

typedef CGAL::Pm_default_point_location<Pm>             Trap_point_location;
typedef CGAL::Pm_naive_point_location<Pm>               Naive_point_location;
typedef CGAL::Pm_walk_along_line_point_location<Pm>     Walk_point_location;
typedef CGAL::Pm_dummy_point_location<Pm>               Dummy_point_location;

inline CGAL::Window_stream & operator<<(CGAL::Window_stream & os, Pmwx & pm)
{
  Pmwx::Edge_iterator ei;
  os << CGAL::BLUE;
  for (ei = pm.edges_begin(); ei != pm.edges_end(); ++ei)
    os << (*ei).curve();
  Pmwx::Vertex_iterator vi;
  os << CGAL::RED;
  for (vi = pm.vertices_begin(); vi != pm.vertices_end(); ++vi)
    os << (*vi).point();
  return os;
}

class Conic_reader
{

public:
  int ReadData(const char * filename, CurveList & curves,
               FormatId format,
               CGAL::Bbox_2 & bbox)
  {
    Curve_2 cv;
    char dummy[256];

    std::ifstream inp(filename);
    if (!inp.is_open()) {
      std::cerr << "Cannot open file " << filename << "!" << std::endl;
      return -1;
    }
    int count;
    inp >> count;
    inp.getline(dummy, sizeof(dummy));
    for (int i = 0; i < count; i++) {
      ReadCurve(inp, cv);
      curves.push_back(cv);
      CGAL::Bbox_2 curve_bbox = cv.bounding_box();
      if (i == 0) bbox = curve_bbox;
      else bbox = bbox + curve_bbox;      
    }
    inp.close();
    return 0;
  }
  
  void ReadCurve(std::ifstream & is, Curve_2 & cv)
  {
      // Read a line from the input file.
      char one_line[128];
      
      skip_comments (is, one_line);
      std::istringstream str_line (one_line);
      
      // Get the arc type.
      char     type;
      bool     is_circle = false;              // Is this a circle.
      Circle_2 circle;
      NT       r, s, t, u, v, w;               // The conic coefficients.
      
      str_line >> type;
      
      // An ellipse (full ellipse or a partial ellipse):
      if (type == 'f' || type == 'F' || type == 'e' || type == 'E')
      {  
          // Read the ellipse (using the format "a b x0 y0"):
          //
          //     x - x0   2      y - y0   2
          //  ( -------- )  + ( -------- )  = 1
          //       a               b
          //
          NT     a, b, x0, y0;
          
          str_line >> a >> b >> x0 >> y0;
          
          NT     a_sq = a*a;
          NT     b_sq = b*b;
          
          if (a == b)
          {
              is_circle = true;
              circle = Circle_2 (Point_2 (x0, y0), a*b, CGAL::CLOCKWISE);
          }
          else
          {
              r = b_sq;
              s = a_sq;
              t = 0;
              u = -2*x0*b_sq;
              v = -2*y0*a_sq;
              w = x0*x0*b_sq + y0*y0*a_sq - a_sq*b_sq;
          }
          
          if (type == 'f' || type == 'F')
          {
              // Create a full ellipse (or circle).
              if (is_circle)
                  cv = Curve_2 (circle);
              else
                  cv = Curve_2 (r, s, t, u, v, w);
              
              return;
          }
      }
      else if (type == 'h' || type == 'H')
      {
          // Read the hyperbola (using the format "a b x0 y0"):
          //
          //     x - x0   2      y - y0   2
          //  ( -------- )  - ( -------- )  = 1
          //       a               b
          //
          NT     a, b, x0, y0;
          
          str_line >> a >> b >> x0 >> y0;
          
          NT     a_sq = a*a;
          NT     b_sq = b*b;
          
          r = b_sq;
          s= -a_sq;
          t = 0;
          u = -2*x0*b_sq;
          v = 2*y0*a_sq;
          w = x0*x0*b_sq - y0*y0*a_sq - a_sq*b_sq;  
      }
      else if (type == 'p' || type == 'P')
      {
          // Read the parabola (using the format "c x0 y0"):
          //
          //                        2
          //  4c*(y - y0) = (x - x0)
          //
          NT     c, x0, y0;
          
          str_line >> c >> x0 >> y0;
          
          r = 1;
          s = 0;
          t = 0;
          u = -2*x0;
          v = -4*c;
          w = x0*x0 + 4*c*y0;
      }
      else if (type == 'c' || type == 'C' || type == 'a' || type == 'A')
      {
          // Read a general conic, given by its coefficients <r,s,t,u,v,w>.
          str_line >> r >> s >> t >> u >> v >> w;
          
          if (type == 'c' || type == 'C')
          {
              // Create a full conic (should work only for ellipses).
              cv = Curve_2 (r, s, t, u, v, w);
              return;
          }
      }
      else if (type == 's' || type == 'S')
      {
          // Read a segment, given by its endpoints (x1,y1) and (x2,y2);
          NT      x1, y1, x2, y2;
          
          str_line >> x1 >> y1 >> x2 >> y2;
          
          Point_2   source (x1, y1);
          Point_2   target (x2, y2);
          Segment_2 segment (source, target);
          
          // Create the segment.
          cv = Curve_2(segment);
          return;
      }
      else if (type == 'i' || type == 'I')
      {
          // Read a general conic, given by its coefficients <r,s,t,u,v,w>.
          str_line >> r >> s >> t >> u >> v >> w;
          
          // Read the approximated source, along with a general conic 
          // <r_1,s_1,t_1,u_1,v_1,w_1> whose intersection with <r,s,t,u,v,w>
          // defines the source.
          NT     r1, s1, t1, u1, v1, w1;
          NT     x1, y1;
          
          str_line >> x1 >> y1;
          str_line >> r1 >> s1 >> t1 >> u1 >> v1 >> w1;
          
          Point_2   app_source (x1, y1);
          
          // Read the approximated target, along with a general conic 
          // <r_2,s_2,t_2,u_2,v_2,w_2> whose intersection with <r,s,t,u,v,w>
          // defines the target.
          NT     r2, s2, t2, u2, v2, w2;
          NT     x2, y2;
          
          str_line >> x2 >> y2;
          str_line >> r2 >> s2 >> t2 >> u2 >> v2 >> w2;
          
          Point_2   app_target (x2, y2);
          
          // Create the conic arc.
          cv = Curve_2 (r, s, t, u, v ,w,
                        app_source, r1, s1, t1, u1, v1, w1,
                        app_target, r2, s2, t2, u2, v2, w2);
          return;
      }
      else
      {
          std::cerr << "Illegal conic type specification: " << type << "."
                    << std::endl;
          return;
      }
      
      // Read the end points of the arc and create it.
      NT    x1, y1, x2, y2;
      
      str_line >> x1 >> y1 >> x2 >> y2;
      
      Point_2 source (x1, y1);
      Point_2 target (x2, y2);
      
      // Create the conic (or circular) arc.
      if (is_circle)
      {
          cv = Curve_2 (circle,
                        source, target);
      }
      else
      {
          cv = Curve_2 (r, s, t, u, v, w,
                        source, target);
      }
      
      return;
  }
    
  void skip_comments( std::ifstream& is, char* one_line )
  {
    while( !is.eof() ){
      is.getline( one_line, 128 );
      if( one_line[0] != '#' ){
	break;
      }
    }  
  }
};

/*! redraw
 */
static void redraw(leda_window * wp, double x0, double y0,
                   double x1, double y1)
{ wp->flush_buffer(x0, y0, x1, y1); }

/*!
 */
int main(int argc, char * argv[])
{
  int verbose = 1;
  CGAL::Bbox_2 bbox;
  FormatId format = FORMAT_INT;
  if (argc != 2) {
    std::cerr << "usage: Conic_arr_from_file filename\n";
    return -1;
  }
  const char * filename = argv[1];
  CurveList curveList;

  // read
  Conic_reader reader;
  int rc = reader.ReadData(filename, curveList, format, bbox);
  if (rc < 0) return rc;
  if (verbose) std::cout << curveList.size() << " curves" << std::endl;
  
  // construct
  Naive_point_location strategy;
  Pmwx pm(&strategy);

  CGAL::Timer t;
  t.start();

  CurveList::const_iterator i;
  for (i = curveList.begin(); i != curveList.end(); i++)
    pm.insert(*i);
 
  t.stop();
  std::cout << "Construction took " 
	    << t.time() << " seconds." << std::endl;
  
  curveList.clear();

  // if map is empty
  if (pm.halfedges_begin() == pm.halfedges_end()) {
      std::cout << std::endl;
      std::cout << "No edges were inserted. Planar map is empty. Exiting.";
      std::cout << std::endl;
      return -1;
  }

  if (verbose) {
    if (!pm.is_valid()) std::cerr << "map invalid!" << std::endl;
    std::cout << "# of vertices: " << pm.number_of_vertices() << std::endl;
    std::cout << "# of halfedges: " << pm.number_of_halfedges() << std::endl;
    std::cout << "# of faces: " << pm.number_of_faces() << std::endl;
  }

  // initialize window
  float x_range = bbox.xmax() - bbox.xmin();
  float y_range = bbox.ymax() - bbox.ymin();
  float width = 640;
  float height = (y_range * width) / x_range;
    
  CGAL::Window_stream * myWindow =
    new CGAL::Window_stream(static_cast<int>(width),
                            static_cast<int>(height),
			    "CGAL - Conic Arcs Arrangement Demo");
  if (!myWindow) return -1;

  float min_range = (x_range < y_range) ? x_range : y_range;
  float x_margin = min_range / 4;
  float y_margin = (height * x_margin) / width;
        
  float x0 = bbox.xmin() - x_margin;
  float x1 = bbox.xmax() + x_margin;
  float y0 = bbox.ymin() - y_margin;
  myWindow->init(x0, x1, y0);   // logical window size 

  myWindow->set_redraw(redraw);
  myWindow->set_mode(CGAL_LEDA_SCOPE::src_mode);
  myWindow->set_node_width(3);
  myWindow->set_point_style(leda_cross_point);
  myWindow->set_line_width(1);
  myWindow->button("finish",10);
  myWindow->display(leda_window::center, leda_window::center);

  // display
  myWindow->set_flush(0);
  (*myWindow) << pm;
  myWindow->set_flush(1);
  myWindow->flush();

  // Point Location Queries
  myWindow->set_status_string("Enter a query point with left mouse button. "
                              "Finish button - exit." );
  (*myWindow) << CGAL::RED;

  Point_2 p;
  Pmwx::Halfedge_handle e;
  
  CGAL::My_Arr_drawer< Pmwx, Pmwx::Ccb_halfedge_circulator,
    Pmwx::Holes_iterator> drawer(*myWindow);
  for (; ; ) {
    double x,y;
    int b = myWindow->read_mouse(x,y);
    if (b == 10) break;
    else
      p = Point_2(x, y);

    (*myWindow) << pm;
    
    Pmwx::Locate_type lt;
    e = pm.locate(p, lt);
      
    //color the face on the screen
    Pmwx::Face_handle f = e->face();
    drawer.draw_face(f);
  }

  delete myWindow;
  return 0;
}

#endif
