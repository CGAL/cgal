#ifndef CONIC_READER_H
#define CONIC_READER_H

#include <CGAL/Bench_parse_args.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <list>
#include <string>

#include "numberType.h"

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
    std::string stringvalues(one_line);
    std::istringstream str_line (stringvalues, std::istringstream::in);
      
    // Get the arc type.
    char type;
    bool is_circle = false;              // Is this a circle.
    CfNT a, b, x0, y0;

#if BENCH_TRAITS != CK_CIRCLE_TRAITS
    CfNT r, s, t, u, v, w;               // The conic coefficients.
#endif

#if BENCH_TRAITS == CONIC_TRAITS || BENCH_TRAITS == CK_CIRCLE_TRAITS || \
    BENCH_TRAITS == CK_CONIC_TRAITS
    Circle_2 circle;
#endif

#if BENCH_TRAITS == CK_CONIC_TRAITS
    Conic_2 conic;
#endif
    
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
      CfRT a_in, b_in, x0_in, y0_in;
      str_line >> a_in >> b_in >> x0_in >> y0_in;
      a = a_in;
      b = b_in;
      x0 = x0_in;
      y0 = y0_in;
          
      if (a == b)
      {
        is_circle = true;
#if BENCH_TRAITS == CONIC_TRAITS || BENCH_TRAITS == CK_CONIC_TRAITS || \
    BENCH_TRAITS == CK_CIRCLE_TRAITS
        circle = Circle_2 (Point_2 (x0, y0), a*b, CGAL::CLOCKWISE);
#endif
      }
      else
      {
#if BENCH_TRAITS != CK_CIRCLE_TRAITS
        CfNT a_sq = a*a;
        CfNT b_sq = b*b;
          
        r = b_sq;
        s = a_sq;
        t = 0;
        u = -2*x0*b_sq;
        v = -2*y0*a_sq;
        w = x0*x0*b_sq + y0*y0*a_sq - a_sq*b_sq;
#if BENCH_TRAITS == CK_CONIC_TRAITS
        conic = Conic_2(r, s, t, u, v, w);
#endif
#endif
      }
          
      if (type == 'f' || type == 'F')
      {
        // Create a full ellipse (or circle).
        if (is_circle) {
#if BENCH_TRAITS == CONIC_TRAITS || BENCH_TRAITS == CK_CIRCLE_TRAITS || \
    BENCH_TRAITS == CK_CONIC_TRAITS
          cv = Curve_2 (circle);
#else
          cv = Curve_2 (x0, y0, a);
#endif              
        } else {
#if BENCH_TRAITS == CK_CIRCLE_TRAITS
          std::cerr << "Skipping Ellipsoid Arc!" << std::endl;
          return false;
#else

#if BENCH_TRAITS == CK_CONIC_TRAITS
          cv = Curve_2 (conic);          
#else
          cv = Curve_2 (r, s, t, u, v, w);
#endif
#endif
        }
        return true;
      }
    }

    else if (type == 'h' || type == 'H')
    {
#if BENCH_TRAITS == CK_CIRCLE_TRAITS || BENCH_TRAITS == CK_CONIC_TRAITS
      std::cerr << "Skipping Hyperbola!" << std::endl;
      return false;
#else
      
      // Read the hyperbola (using the format "a b x0 y0"):
      //
      //     x - x0   2      y - y0   2
      //  ( -------- )  - ( -------- )  = 1
      //       a               b
      //
      CfNT a, b, x0, y0;
      CfRT a_in, b_in, x0_in, y0_in;
      str_line >> a_in >> b_in >> x0_in >> y0_in;
      a = a_in;
      b = b_in;
      x0 = x0_in;
      y0 = y0_in;
          
      CfNT a_sq = a*a;
      CfNT b_sq = b*b;
          
      r = b_sq;
      s= -a_sq;
      t = 0;
      u = -2*x0*b_sq;
      v = 2*y0*a_sq;
      w = x0*x0*b_sq - y0*y0*a_sq - a_sq*b_sq;
#endif      
    }
    else if (type == 'p' || type == 'P')
    {
#if BENCH_TRAITS == CK_CIRCLE_TRAITS || BENCH_TRAITS == CK_CONIC_TRAITS
      std::cerr << "Skipping Parabola!" << std::endl;
      return false;
#else
      
      // Read the parabola (using the format "c x0 y0"):
      //
      //                        2
      //  4c*(y - y0) = (x - x0)
      //
      CfNT c, x0, y0;
      CfRT c_in, x0_in, y0_in;
      str_line >> c_in >> x0_in >> y0_in;
      c = c_in;
      x0 = x0_in;
      y0 = y0_in;
          
      r = 1;
      s = 0;
      t = 0;
      u = -2*x0;
      v = -4*c;
      w = x0*x0 + 4*c*y0;
#endif
    }
    else if (type == 'c' || type == 'C' || type == 'a' || type == 'A')
    {
#if BENCH_TRAITS == CK_CIRCLE_TRAITS
      std::cerr << "Skipping Generl Conic!" << std::endl;
      return false;
#else

      // Read a general conic, given by its coefficients <r,s,t,u,v,w>.
      CfRT r_in, s_in, t_in, u_in, v_in, w_in;
      str_line >> r_in >> s_in >> t_in >> u_in >> v_in >> w_in;
      r = r_in;
      s = s_in;
      t = t_in;
      u = u_in;
      v = v_in;
      w = w_in;
          
      if (type == 'c' || type == 'C')
      {
        // Create a full conic (should work only for ellipses).
#if BENCH_TRAITS == CK_CONIC_TRAITS
        conic = Conic_2(r, s, t, u, v, w);
        if (conic.is_ellipse())
          cv = Curve_2 (conic);
        else {
          std::cerr << "Skipping Non Ellipse" << std::endl;
          return false;
        }
#else
        cv = Curve_2 (r, s, t, u, v, w);
#endif
        return true;
      }
#endif
    }
    else if (type == 's' || type == 'S')
    {
#if BENCH_TRAITS == CK_CIRCLE_TRAITS || BENCH_TRAITS == CK_CONIC_TRAITS
      std::cerr << "Skipping Segment!" << std::endl;
      return false;
#else
      
      // Read a segment, given by its endpoints (x1,y1) and (x2,y2);
      CfNT x1, y1, x2, y2;
      CfRT x1_in, y1_in, x2_in, y2_in;
      str_line >> x1_in >> y1_in >> x2_in >> y2_in;
      x1 = x1_in;
      y1 = y1_in;
      x2 = x2_in;
      y2 = y2_in;
      
      // Create the segment.
#if BENCH_TRAITS == CORE_CONIC_TRAITS
      cv = Curve_2(x1, y1, x2, y2);
#else
      Point_2 source (x1, y1);
      Point_2 target (x2, y2);
      Segment_2 segment (source, target);
      cv = Curve_2(segment);
#endif
      return true;
#endif
    }
    else if (type == 'i' || type == 'I')
    {
#if BENCH_TRAITS == CK_CIRCLE_TRAITS || BENCH_TRAITS == CK_CONIC_TRAITS
      std::cerr << "Skipping Approximate Arc!" << std::endl;
      return false;
#else
      
      // Read a general conic, given by its coefficients <r,s,t,u,v,w>.
      CfRT r_in, s_in, t_in, u_in, v_in, w_in;
      str_line >> r_in >> s_in >> t_in >> u_in >> v_in >> w_in;
      r = r_in;
      s = s_in;
      t = t_in;
      u = u_in;
      v = v_in;
      w = w_in;
          
      // Read the approximated source, along with a general conic 
      // <r_1,s_1,t_1,u_1,v_1,w_1> whose intersection with <r,s,t,u,v,w>
      // defines the source.
      CfNT r1, s1, t1, u1, v1, w1;
      CfNT x1, y1;
          
      CfRT x1_in, y1_in;
      str_line >> x1_in >> y1_in;
      x1 = x1_in;
      y1 = y1_in;
      
      CfRT r1_in, s1_in, t1_in, u1_in, v1_in, w1_in;
      str_line >> r1_in >> s1_in >> t1_in >> u1_in >> v1_in >> w1_in;
      r1 = r1_in;
      s1 = s1_in;
      t1 = t1_in;
      u1 = u1_in;
      v1 = v1_in;
      w1 = w1_in;
      
      // Read the approximated target, along with a general conic 
      // <r_2,s_2,t_2,u_2,v_2,w_2> whose intersection with <r,s,t,u,v,w>
      // defines the target.
      CfNT r2, s2, t2, u2, v2, w2;
      CfNT x2, y2;
      
      CfRT x2_in, y2_in;
      str_line >> x2_in >> y2_in;
      x2 = x2_in;
      y2 = y2_in;

      CfRT r2_in, s2_in, t2_in, u2_in, v2_in, w2_in;
      str_line >> r2_in >> s2_in >> t2_in >> u2_in >> v2_in >> w2_in;
      r2 = r2_in;
      s2 = s2_in;
      t2 = t2_in;
      u2 = u2_in;
      v2 = v2_in;
      w2 = w2_in;
      
      // Create the conic arc.
#if BENCH_TRAITS == CORE_CONIC_TRAITS
      std::cerr << "Not implemented!" << std::endl;
      return false;
#else
      Point_2 app_source (x1, y1);
      Point_2 app_target (x2, y2);
      cv = Curve_2 (r, s, t, u, v ,w,
                    app_source, r1, s1, t1, u1, v1, w1,
                    app_target, r2, s2, t2, u2, v2, w2);
#endif
      return true;
#endif
    }
    else
    {
      std::cerr << "Illegal conic type specification: " << type << "."
                << std::endl;
      return false;
    }

    // Read the end points of the arc and create it.
    CfNT x1, y1, x2, y2;
    CfRT x1_in, y1_in, x2_in, y2_in;
    str_line >> x1_in >> y1_in >> x2_in >> y2_in;
    x1 = x1_in;
    y1 = y1_in;
    x2 = x2_in;
    y2 = y2_in;
      
    Point_2 source (x1, y1);
    Point_2 target (x2, y2);
      
    // Create the conic (or circular) arc.
    if (is_circle)
    {
#if BENCH_TRAITS != CK_CIRCLE_TRAITS
#if BENCH_TRAITS == CORE_CONIC_TRAITS
      cv = Curve_2 (x0, y0, a, CGAL::CLOCKWISE, source, target);
#else
      cv = Curve_2 (circle, source, target);
#endif
    }
    else
    {
#if BENCH_TRAITS == CK_CONIC_TRAITS
      if (conic.is_ellipse()) {
        cv = Curve_2 (conic, source, target);
      } else {
        std::cerr << "Skipping Non Ellipse" << std::endl;
        return false;
      }
#elif BENCH_TRAITS == CORE_CONIC_TRAITS
      cv = Curve_2 (r, s, t, u, v, w, CGAL::CLOCKWISE, source, target);
#else
      cv = Curve_2 (r, s, t, u, v, w, source, target);
#endif
#endif
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
