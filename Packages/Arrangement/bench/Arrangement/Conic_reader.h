#ifndef CONIC_READER_H
#define CONIC_READER_H

#include <CGAL/leda_real.h>
#include <CGAL/Bench_parse_args.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <list>
#include <string>

#include "numberType.h"

template <class Traits>
class Conic_reader
{

public:
  typedef typename Traits::Point_2    Point_2;
  typedef typename Traits::Curve_2    Curve_2;
  typedef typename Traits::Circle_2   Circle_2;
  typedef typename Traits::Segment_2  Segment_2;
  typedef typename Traits::X_curve_2  X_curve_2;
  typedef typename std::list<Curve_2> CurveList;

  int ReadData(const char * filename, CurveList & curves,
               CGAL::Bench_parse_args::FormatId format,
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
      std::string stringvalues(one_line);
      std::istringstream str_line (stringvalues, std::istringstream::in);
      
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

#endif
