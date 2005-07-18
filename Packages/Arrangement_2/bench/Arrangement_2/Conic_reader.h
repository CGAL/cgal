#ifndef CONIC_READER_H
#define CONIC_READER_H

#include <CGAL/Bench_parse_args.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <list>
#include <string>

#include "number_type.h"

#if BENCH_TRAITS == EXACUS_CONIC_TRAITS
#include <LiS/file_io.h>
#endif

#if BENCH_NT == CORE_EXPR_NT
#include <CORE/BigInt.h>
#endif

template <class Traits>
class Conic_reader
{
public:
  typedef typename Traits::Curve_2              Curve_2;
  typedef typename Traits::X_monotone_curve_2   X_monotone_curve_2;

#if BENCH_TRAITS == CK_CIRCLE_TRAITS || BENCH_TRAITS == CK_CONIC_TRAITS
  typedef typename Traits::Kernel               Curved_k;
  typedef typename Curved_k::Linear_kernel      Linear_k;
  typedef typename Curved_k::Circle_2           Circle_2;
  typedef typename Curved_k::Point_2            Point_2;
#if BENCH_NT == LAZY_CGAL_GMPQ_NT
  typedef CGAL::Gmpq                            CfNT;
#else
  typedef typename Linear_k::FT                 CfNT;
#endif
#if 0
  typedef typename Linear_k::RT                 CfRT;
#elif BENCH_NT == CORE_EXPR_NT
  typedef CORE::BigInt                          CfRT;
#else
  typedef RT                                    CfRT;
#endif
  
#if BENCH_TRAITS == CK_CONIC_TRAITS
  typedef typename Curved_k::Conic_2            Conic_2;  
#endif
  
#else
  typedef typename Traits::Point_2              Point_2;

#if BENCH_TRAITS == CORE_CONIC_TRAITS
  typedef typename Traits::CfNT                 CfNT;
  typedef typename Traits::CfNT                 CfRT;
  typedef typename Traits::CoNT                 CoNT;
  typedef typename Traits::Int_point_2          Int_point_2;
  typedef typename Traits::Int_segment_2        Int_segment_2;
  typedef typename Traits::Int_circle_2         Int_circle_2;
  typedef typename Traits::Int_line_2           Int_line_2;

#elif BENCH_TRAITS != EXACUS_CONIC_TRAITS
  typedef typename Traits::Circle_2             Circle_2;
  typedef typename Traits::Segment_2            Segment_2;
  typedef typename Traits::NT                   CfNT;
#if 0
  typedef typename Traits::RT                   CfRT;
#else
  typedef typename Traits::NT                   CfRT;
#endif
#endif
#endif
  
  template<class OutputIterator>
  int read_data(const char * filename, OutputIterator curves_out,
                CGAL::Bench_parse_args::FormatId format, CGAL::Bbox_2 & bbox)
  {
#if BENCH_TRAITS == EXACUS_CONIC_TRAITS
    Curve_2 cv;
      
    std::list< Curve_2 > curves;
    bool success = LiS::read_file(filename, curves);
    if (!success) {
      return 0;
    }
    std::copy(curves.begin(), curves.end(), curves_out);
#if 0
    // TODO set these boxes
    for (typename std::list< Curve_2 >::iterator it = curves.begin();
         it != curves.end();
         it++)
    {
      CGAL::Bbox_2 curve_bbox = cv.bbox();
      if (i == 0) bbox = curve_bbox;
      else bbox = bbox + curve_bbox;
    }
#else
    bbox = CGAL::Bbox_2(-1000, -1000, 1000, 1000);
#endif
    return 0;
      
#else
      
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
#endif
  }

#if BENCH_TRAITS != EXACUS_CONIC_TRAITS
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
#endif
};

#endif
