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
#include <CGAL/Arr_segment_traits_2.h>

#include <CGAL/Pm_default_dcel.h>
#include <CGAL/Planar_map_2.h>
#include <CGAL/Pm_with_intersections.h>
#include <CGAL/IO/Window_stream.h>
#include <CGAL/IO/Pm_iostream.h>
#include <CGAL/IO/Pm_Window_stream.h>

#include <CGAL/Pm_trapezoid_ric_point_location.h>
#include <CGAL/Pm_walk_along_line_point_location.h>
#include <CGAL/Pm_naive_point_location.h>
#include <CGAL/Pm_dummy_point_location.h>

#include <stdlib.h>
#include <iostream>
#include <fstream>
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

typedef CGAL::Quotient<CGAL::MP_Float>                  NT;
typedef CGAL::Cartesian<NT>                             Kernel;
typedef CGAL::Arr_segment_traits_2<Kernel>              Traits;
typedef CGAL::Pm_default_dcel<Traits>                   Dcel;
typedef CGAL::Planar_map_2<Dcel,Traits>                 Pm;
typedef CGAL::Planar_map_with_intersections_2<Pm>       Pmwx;
typedef Traits::Point_2                                 Point;
typedef Traits::X_monotone_curve_2                      Curve;
typedef std::list<Curve>                                CurveList;

typedef CGAL::Pm_trapezoid_ric_point_location<Pm>             Trap_point_location;
typedef CGAL::Pm_naive_point_location<Pm>               Naive_point_location;
typedef CGAL::Pm_walk_along_line_point_location<Pm>     Walk_point_location;
typedef CGAL::Pm_dummy_point_location<Pm>               Dummy_point_location;

#if defined(LEDA_KERNEL)
inline CGAL::Window_stream & operator<<(CGAL::Window_stream & os,
                                        const Point & p)
{
  os << leda_point(p.xcoordD(), p.ycoordD());
  return os;
}

inline CGAL::Window_stream & operator<<(CGAL::Window_stream & os,
                                        const Kernel::Segment_2 & seg)
{
  os << leda_segment(seg.xcoord1D(), seg.ycoord1D(),
                     seg.xcoord2D(), seg.ycoord2D());
  return os;
}
#endif

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

/*!
 */
class Segment_reader
{
public:
  int ReadData(const char * filename, CurveList & curveList,
               FormatId format, CGAL::Bbox_2 & bbox)
  {
    std::ifstream inp(filename);
    if (!inp.is_open()) {
      std::cerr << "Cannot open file " << filename << "!" << std::endl;
      return -1;
    }
    int count;
    inp >> count;
    
    int i;
    for (i = 0; i < count; i++) {
      NT x0, y0, x1, y1;
      if (format == FORMAT_RAT) {
        inp >> x0 >> y0 >> x1 >> y1;
      } else if (format == FORMAT_INT) {
        int ix0, iy0, ix1, iy1;
        inp >> ix0 >> iy0 >> ix1 >> iy1;
        x0 = ix0; y0 = iy0; x1 = ix1; y1 = iy1;
      } else if (format == FORMAT_FLT) {
        float ix0, iy0, ix1, iy1;
        inp >> ix0 >> iy0 >> ix1 >> iy1;
        x0 = ix0; y0 = iy0; x1 = ix1; y1 = iy1;
      } else {
        std::cerr << "Illegal format!" << std::endl;
        return -1;
      }

      Point p1(x0, y0);
      Point p2(x1, y1);

      // if (p1 == p2) continue;
      Curve curve(p1, p2);
      curveList.push_back(curve);

      // Update the bounding box of the arrangement.
#if defined(USE_LEDA)
      double xmin, ymin, xmax, ymax;
      if (p1.xcoord() < p2.xcoord()) {
        xmin = CGAL::to_double(p1.xcoord());
        xmax = CGAL::to_double(p2.xcoord());
      } else {
        xmin = CGAL::to_double(p2.xcoord());
        xmax = CGAL::to_double(p1.xcoord());
      }
      if (p1.ycoord() < p2.ycoord()) {
        ymin = CGAL::to_double(p1.ycoord());
        ymax = CGAL::to_double(p2.ycoord());
      } else {
        ymin = CGAL::to_double(p2.ycoord());
        ymax = CGAL::to_double(p1.ycoord());
      }
      
      CGAL::Bbox_2 curve_bbox(xmin, ymin, xmax, ymax);
#else
      CGAL::Bbox_2 curve_bbox = curve.bbox();
#endif
      if (i == 0) bbox = curve_bbox;
      else bbox = bbox + curve_bbox;
    }
    inp.close();
    return 0;
  }
};

/*! redraw
 */
static void redraw(CGAL::Window_stream * wp, double x0, double y0,
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
    std::cerr << "usage: Segment_arr_from_file filename\n";
    return -1;
  }
  const char * filename = argv[1];
  CurveList curveList;

  // read
  Segment_reader reader;
  int rc = reader.ReadData(filename, curveList, format, bbox);
  if (rc < 0) return rc;
  if (verbose) std::cout << curveList.size() << " curves" << std::endl;
  
  // construct
  Naive_point_location strategy;
  Pmwx pm(&strategy);
  pm.insert(curveList.begin(), curveList.end());
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
                            static_cast<int>(height));
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

  // (*myWindow) << pm;

  // Point Location Queries
  myWindow->set_status_string("Enter a query point with left mouse button. "
                              "Finish button - exit." );
  (*myWindow) << CGAL::RED;

  Point p;
  Pmwx::Halfedge_handle e;
  
  CGAL::My_Arr_drawer< Pmwx, Pmwx::Ccb_halfedge_circulator,
    Pmwx::Holes_iterator> drawer(*myWindow);
  for (; ; ) {
    double x,y;
    int b = myWindow->read_mouse(x,y);
    if (b == 10) break;
    else
      p = Point(x, y);

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
