#ifndef CONIC_READER_H
#define CONIC_READER_H

#include <CGAL/Bench_parse_args.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <list>
#include <string>

#include "number_type.h"

template <class Traits>
class Conic_reader
{
public:
  typedef typename Traits::Curve_2                      Curve_2;
  typedef typename Traits::X_monotone_curve_2           X_monotone_curve_2;
  typedef typename Traits::Point_2                      Point_2;

  typedef typename Traits::Rat_kernel                   Rat_kernel;
  typedef typename Rat_kernel::FT                       Rational;
  
  typedef typename Rat_kernel::Point_2                  Rat_point_2;
  typedef typename Rat_kernel::Segment_2                Rat_segment_2;
  typedef typename Rat_kernel::Line_2                   Rat_line_2;
  typedef typename Rat_kernel::Circle_2                 Rat_circle_2;
  
  template<class OutputIterator>
  int read_data(const char * filename, OutputIterator curves_out,
                CGAL::Bench_parse_args::FormatId format, CGAL::Bbox_2 & bbox)
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
      if (read_curve(inp, cv)) {
        ++curves_out = cv;
        CGAL::Bbox_2 curve_bbox = cv.bbox();
        if (i == 0) bbox = curve_bbox;
        else bbox = bbox + curve_bbox;
      }
    }
    inp.close();
    return 0;
  }

  /*! */
  bool read_curve(std::ifstream & is, Curve_2 & cv)
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
      CORE::BigInt x1, y1, x2, y2;
    
      str_line >> x1 >> y1 >> x2 >> y2;
    
      Rat_point_2   p1(x1, y1), p2(x2, y2);
      Rat_segment_2 seg (p1, p2);
      cv = Curve_2 (seg);
    }
    else if (type == 'c' || type == 'C')
    {
      // Construct a full circle. The line should have the format:
      //   c <x0> <y0> <R_sq>
      // where (x0, y0) is the center of the circle and R_sq is its squared
      // radius.
      CORE::BigInt    x0, y0, R_sq;
      
      str_line >> x0 >> y0 >> R_sq;
      Rat_point_2   p0(x0, y0);
      Rat_circle_2  circ(p0, R_sq);
      cv = Curve_2 (circ);
    }
    else if (type == 't' || type == 'T')
    {
      // Construct a circular arc. The line should have the format:
      //   t <x1> <y1> <x2> <y2> <x3> <y3>
      // where (x1, y1), (x2, y2) and (x3, y3) define the arc.
      CORE::BigInt    x1, y1, x2, y2, x3, y3;
      
      str_line >> x1 >> y1 >> x2 >> y2 >> x3 >> y3;
      Rat_point_2   p1(x1, y1), p2(x2, y2), p3(x3, y3);
      cv = Curve_2 (p1, p2, p3);
    }
    else if (type == 'f' || type == 'F')
    {
      // Construct a full conic curve. The line should have the format:
      //   c <r> <s> <t> <u> <v> <w>
      // where r, s, t, u, v, w define the conic equation.
      CORE::BigInt    r, s, t, u, v, w;
      
      str_line >> r >> s >> t >> u >> v >> w;
      cv = Curve_2 (r, s, t, u, v, w);
    }
    else if (type == 'a' || type == 'A')
    {
      // Construct a conic arc. The line should have the format:
      //   c <r> <s> <t> <u> <v> <w> <orient> <x1> <y1> <x2> <y2>
      // where r, s, t, u, v, w define the conic equation, while (x1, y1)
      // and (x2, y2) are the arc's endpoints.
      CORE::BigInt    r, s, t, u, v, w;
      
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
      Rational    x1, y1, x2, y2;
      
      str_line >> num;
      x1 = Rational(num);
      str_line >> num;
      y1 = Rational(num);
      
      str_line >> num;
      x2 = Rational(num);
      str_line >> num;
      y2 = Rational(num);
      
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
      CORE::BigInt    r, s, t, u, v, w;
      CORE::BigInt    a, b, c;
      
      str_line >> r >> s >> t >> u >> v >> w >> a >> b >> c;
      Rat_line_2    line (a, b, c);
      cv = Curve_2 (r, s, t, u, v, w, line);
    }
    else if (type == 'q' || type == 'Q')
    {
      // Construct a circular arc. The line should have the format:
      //   t <x1> <y1> <x2> <y2> <x3> <y3> <x4> <y4> <x5> <y5>
      // where (x1, y1), (x2, y2), (x3, y3), (x4, y4) and (x5, y5) define the 
      // arc.
      CORE::BigInt    x1, y1, x2, y2, x3, y3, x4, y4, x5, y5;
      
      str_line >> x1 >> y1 >> x2 >> y2 >> x3 >> y3 >> x4 >> y4 >> x5 >> y5;
      Rat_point_2   p1(x1, y1), p2(x2, y2), p3(x3, y3), p4(x4, y4), p5(x5, y5);
      cv = Curve_2 (p1, p2, p3, p4, p5);
    }
    else
    {
      std::cerr << "Illegal conic type specification: " << type << "."
                << std::endl;
      return false;
    }
    
    return true;
  }

  /*! */
  void skip_comments( std::ifstream& is, char* one_line )
  {
    while (!is.eof()) {
      is.getline(one_line, 128);
      if (one_line[0] != '#') break;
    }
  }
};

#endif
