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
#include <CORE/BigInt.h>
#include <CGAL/CORE_Expr.h>
#include <CGAL/Arr_conic_traits_2.h>
#include <CGAL/Timer.h>

#include <CGAL/Pm_default_dcel.h>
#include <CGAL/Planar_map_2.h>
#include <CGAL/Pm_with_intersections.h>
#include <CGAL/IO/Window_stream.h>
#include <CGAL/IO/Pm_iostream.h>
#include <CGAL/IO/Pm_Window_stream.h>
#include <CGAL/IO/Conic_arc_2_Window_stream.h>

#include <CGAL/Pm_trapezoid_ric_point_location.h>
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

typedef CORE::BigInt                                CfNT;
typedef CGAL::Cartesian<CfNT>                       Int_kernel;
typedef Int_kernel::Point_2                         Int_point_2;
typedef Int_kernel::Segment_2                       Int_segment_2;
typedef Int_kernel::Line_2                          Int_line_2;
typedef Int_kernel::Circle_2                        Int_circle_2;

typedef CORE::Expr                                      CoNT;
typedef CGAL::Cartesian<CoNT>                           Alg_kernel;

typedef CGAL::Arr_conic_traits_2<Int_kernel,Alg_kernel> Traits_2;
typedef CGAL::Pm_default_dcel<Traits_2>                 Dcel;
typedef CGAL::Planar_map_2<Dcel,Traits_2>               Pm_2;
typedef CGAL::Planar_map_with_intersections_2<Pm_2>     Pmwx_2;

typedef Traits_2::Point_2                               Point_2;
typedef Traits_2::Curve_2                               Curve_2;
typedef Traits_2::X_monotone_curve_2                    X_monotone_curve_2;
typedef std::list<Curve_2>                              CurveList;

typedef CGAL::Pm_trapezoid_ric_point_location<Pm_2>     Trap_point_location;
typedef CGAL::Pm_naive_point_location<Pm_2>             Naive_point_location;
typedef CGAL::Pm_walk_along_line_point_location<Pm_2>   Walk_point_location;
typedef CGAL::Pm_dummy_point_location<Pm_2>             Dummy_point_location;

/*!
 * Draw the arrangement.
 */
inline CGAL::Window_stream & operator<<(CGAL::Window_stream & os, Pmwx_2 & pm)
{
  // Draw all edges.
  Pmwx_2::Edge_iterator ei;

  os << CGAL::BLUE;
  for (ei = pm.edges_begin(); ei != pm.edges_end(); ++ei)
    os << (*ei).curve();

  // Draw all vertices.
  Pmwx_2::Vertex_iterator vi;
  os << CGAL::RED;
  for (vi = pm.vertices_begin(); vi != pm.vertices_end(); ++vi)
    os << (*vi).point();

  return (os);
}

/*!
 * A class for reading conic arcs from an input file.
 */
class Conic_reader
{

public:
  
  /*!
   * Read a list of curves from a file.
   */
  int ReadData(const char * filename, CurveList & curves,
               FormatId format,
               CGAL::Bbox_2 & bbox)
  {
    Curve_2 cv;
    char dummy[256];

    // Open the input file.
    std::ifstream inp(filename);
    
    if (!inp.is_open()) 
    {
      std::cerr << "Cannot open the input file <" 
		<< filename << ">." << std::endl;
      return (-1);
    }

    // Read the number of curves.
    int count;
    inp >> count;
    inp.getline(dummy, sizeof(dummy));

    // Read the curves.
    for (int i = 0; i < count; i++) 
    {
      ReadCurve(inp, cv);
      curves.push_back(cv);
      CGAL::Bbox_2 curve_bbox = cv.bbox();
      
      if (i == 0) 
	bbox = curve_bbox;
      else 
	bbox = bbox + curve_bbox;      
    }
    inp.close();
    return 0;
  }

  /*!
   * Read the next curve from the file.
   */  
  void ReadCurve(std::ifstream & is, Curve_2 & cv)
  {
    // Read a line from the input file.
    char one_line[128];
      
    skip_comments (is, one_line);
    std::istringstream str_line (one_line);
      
    // Read the arc type and act accordingly.
    char     type;
      
    str_line >> type;
      
    if (type == 's' || type == 'S')
    {
      // Construct a line segment. The line should have the format:
      //   s <x1> <y1> <x2> <y2>
      // where (x1, y1), (x2, y2) are the endpoints of a segment.
      CfNT    x1, y1, x2, y2;

      str_line >> x1 >> y1 >> x2 >> y2;

      Int_point_2   p1(x1, y1), p2(x2, y2);
      Int_segment_2 seg (p1, p2);

      cv = Curve_2 (seg);
    }
    else if (type == 'c' || type == 'C')
    {
      // Construct a full circle. The line should have the format:
      //   c <x0> <y0> <R_sq>
      // where (x0, y0) is the center of the circle and R_sq is its squared
      // radius.
      CfNT    x0, y0, R_sq;

      str_line >> x0 >> y0 >> R_sq;

      Int_point_2   p0(x0, y0);
      Int_circle_2  circ(p0, R_sq);

      cv = Curve_2 (circ);
    }
    else if (type == 't' || type == 'T')
    {
      // Construct a circular arc. The line should have the format:
      //   t <x1> <y1> <x2> <y2> <x3> <y3>
      // where (x1, y1), (x2, y2) and (x3, y3) define the arc.
      CfNT    x1, y1, x2, y2, x3, y3;

      str_line >> x1 >> y1 >> x2 >> y2 >> x3 >> y3;

      Int_point_2   p1(x1, y1), p2(x2, y2), p3(x3, y3);

      cv = Curve_2 (p1, p2, p3);
    }
    else if (type == 'f' || type == 'F')
    {
      // Construct a full conic curve. The line should have the format:
      //   c <r> <s> <t> <u> <v> <w>
      // where r, s, t, u, v, w define the conic equation.
      CfNT    r, s, t, u, v, w;

      str_line >> r >> s >> t >> u >> v >> w;

      cv = Curve_2 (r, s, t, u, v, w);
    }
    else if (type == 'a' || type == 'A')
    {
      // Construct a conic arc. The line should have the format:
      //   c <r> <s> <t> <u> <v> <w> <orient> <x1> <y1> <x2> <y2>
      // where r, s, t, u, v, w define the conic equation, while (x1, y1)
      // and (x2, y2) are the arc's endpoints.
      CfNT    r, s, t, u, v, w;

      str_line >> r >> s >> t >> u >> v >> w;

      // Read the orientation.
      int               i_orient;
      CGAL::Orientation orient;
      
      str_line >> i_orient;
      if (i_orient > 0)
	orient = CGAL::COUNTERCLOCKWISE;
      else if (i_orient < 0)
	orient = CGAL::CLOCKWISE;
      else
	orient = CGAL::COLLINEAR;

      // Read the end points of the arc and create it.
      // Notice we read the coordinates as strings, then we convert them to 
      // the CoNT type, as we do not want to initialize CoNT from a double.
      char    num[50];
      CoNT    x1, y1, x2, y2;
      
      str_line >> num;
      x1 = CoNT(num);
      str_line >> num;
      y1 = CoNT(num);

      str_line >> num;
      x2 = CoNT(num);
      str_line >> num;
      y2 = CoNT(num);
      
      Point_2 ps (x1, y1);
      Point_2 pt (x2, y2);

      cv = Curve_2 (r, s, t, u, v, w, orient, ps ,pt);
    }
    else if (type == 'l' || type == 'L')
    {
      // Construct a conic arc. The line should have the format:
      //   c <r> <s> <t> <u> <v> <w> <a> <b> <c>
      // where r, s, t, u, v, w define the conic equation and a, b, c define
      // a line that intersects it.
      CfNT    r, s, t, u, v, w;
      CfNT    a, b, c;

      str_line >> r >> s >> t >> u >> v >> w >> a >> b >> c;

      Int_line_2    line (a, b, c);

      cv = Curve_2 (r, s, t, u, v, w, line);
    }
    else if (type == 'q' || type == 'Q')
    {
      // Construct a circular arc. The line should have the format:
      //   t <x1> <y1> <x2> <y2> <x3> <y3> <x4> <y4> <x5> <y5>
      // where (x1, y1), (x2, y2), (x3, y3), (x4, y4) and (x5, y5) define the 
      // arc.
      CfNT    x1, y1, x2, y2, x3, y3, x4, y4, x5, y5;

      str_line >> x1 >> y1 >> x2 >> y2 >> x3 >> y3 >> x4 >> y4 >> x5 >> y5;

      Int_point_2   p1(x1, y1), p2(x2, y2), p3(x3, y3), p4(x4, y4), p5(x5, y5);

      cv = Curve_2 (p1, p2, p3, p4, p5);
    }
    else
    {
      std::cerr << "Illegal conic type specification: " << type << "."
		<< std::endl;
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

/*! 
 * Redraw
 */
static void redraw(leda_window * wp, double x0, double y0,
                   double x1, double y1)
{ 
  wp->flush_buffer(x0, y0, x1, y1);
}

/*!
 * The main program.
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

  // Read the conic arcs.
  Conic_reader reader;
  int          rc = reader.ReadData (filename, curveList, format, bbox);

  if (rc < 0) 
    return (rc);

  if (verbose)
    std::cout << "Read " << curveList.size() << " curves." << std::endl;
  
  // Construct the arrangement.
  Naive_point_location      strategy;
  Pmwx_2                    pm (&strategy);
  CurveList::const_iterator it;

  CGAL::Timer t;
  t.start();

  for (it = curveList.begin(); it != curveList.end(); it++)
    pm.insert(*it);
 
  t.stop();
  std::cout << "Construction took " 
	    << t.time() << " seconds." << std::endl;
  
  curveList.clear();

  // Stop here if the map is empty.
  if (pm.halfedges_begin() == pm.halfedges_end()) 
  {
      std::cout << std::endl;
      std::cout << "No edges were inserted. Planar map is empty. Exiting.";
      std::cout << std::endl;
      return (-1);
  }

  if (verbose) 
  {
    if (!pm.is_valid()) std::cerr << "map invalid!" << std::endl;
    std::cout << "# of vertices: " << pm.number_of_vertices() << std::endl;
    std::cout << "# of halfedges: " << pm.number_of_halfedges() << std::endl;
    std::cout << "# of faces: " << pm.number_of_faces() << std::endl;
  }

  // Initialize the window.
  float x_range = bbox.xmax() - bbox.xmin();
  float y_range = bbox.ymax() - bbox.ymin();
  float width = 640;
  float height = (y_range * width) / x_range;
    
  CGAL::Window_stream * myWindow =
    new CGAL::Window_stream(static_cast<int>(width),
                            static_cast<int>(height),
			    "CGAL - Conic Arcs Arrangement Demo");
  if (!myWindow) 
    return (-1);

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

  // Display the arrangement.
  myWindow->set_flush(0);
  (*myWindow) << pm;
  myWindow->set_flush(1);
  myWindow->flush();

  // Answer point-location queries.
  myWindow->set_status_string("Enter a query point with left mouse button. "
                              "Finish button - exit." );
  (*myWindow) << CGAL::RED;

  Point_2 p;
  Pmwx_2::Halfedge_handle e;
  
  CGAL::My_Arr_drawer< Pmwx_2, Pmwx_2::Ccb_halfedge_circulator,
    Pmwx_2::Holes_iterator> drawer(*myWindow);
  
  while (true)
  {
    double x,y;
    int b = myWindow->read_mouse(x,y);
    if (b == 10) break;
    else
      p = Point_2(x, y);

    (*myWindow) << pm;
    
    Pmwx_2::Locate_type lt;
    e = pm.locate(p, lt);
      
    // Color the face on the screen.
    Pmwx_2::Face_handle f = e->face();
    drawer.draw_face(f);
  }

  delete myWindow;
  return 0;
}

#endif
