#ifndef CONIC_READER_HPP
#define CONIC_READER_HPP

#include <CGAL/basic.h>
#include <CGAL/Bbox_2.h>
#include <CGAL/Benchmark/config.hpp>
#include <CGAL/Benchmark/Benchmark_visitor.hpp>
#include <CGAL/Benchmark/benchmark_format.hpp>
#include <CGAL/number_utils_classes.h>

#include <iostream>
#include <string>

#include "number_type.hpp"
#include "kernel_type.hpp"
#include "Point_parser_visitor.hpp"

namespace cb = CGAL::benchmark;

template <>
CORE::BigInt lexical_cast<CORE::BigInt>(std::string & str)
{
  std::stringstream ss(str);
  CORE::BigInt result;
  ss >> result;
  return result;
}

/*! A reader of conic curves and arcs of conic curves */ 
template <class Traits>
class Conic_reader {
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
  
  /*! A parser of bff files that contains conic curves or arcs of conic
   * curves
   */
  template <typename OutputIterator>  
  class Conic_parser_visitor :
    public Point_parser_visitor<Rat_kernel, Rat_point_2, CORE::BigInt> {
  private:
    /*! The iterator of the output container */
    OutputIterator & m_output_iterator;

    /*! Big ints used to store coefficients */
    std::list<CORE::BigInt> m_bigints;
    
  public:
    /*! Constructor */
    Conic_parser_visitor(OutputIterator & oi, CGAL::Bbox_2 & bbox) :
      Point_parser_visitor<Rat_kernel, Rat_point_2, CORE::BigInt>(bbox),
      m_output_iterator(oi)
    {}
    
    /*! Accept only unbounded lines for Arrangements */ 
    virtual void accept_classification(std::string problem,
                                       std::string geom, 
                                       std::string clas,
                                       std::string family,
                                       std::string instance,
                                       std::string release)
    {
      if (problem != "Arrangement")
        this->error_handler( "classification error");
      if (geom != "Conics")
        this->error_handler( "classification error" );
      if (clas != "BoundedArcs")
        this->error_handler( "classification error" );
    }

    /*! Process start line */
    virtual void begin_line_segment_2() { this->m_points.clear(); }

    /*! Process end line
     * Add a curve to the container the output iterator refers to
     */
    virtual void end_line_segment_2()
    {
      Rat_segment_2 seg (this->m_points.front(), this->m_points.back());
      Curve_2 conic(seg);
      ++m_output_iterator = conic;
    }

    /*! Process start circle */
    virtual void begin_circle_2() { this->m_points.clear(); }

    virtual void accept_integer(std::string s)
    {
      CORE::BigInt integer = lexical_cast<CORE::BigInt>(s);
      m_bigints.push_back(integer);
    }
    
    /*! Process end conic
     * Add a curve to the container the output iterator refers to
     */
    virtual void end_circle_2()
    {
      std::list<typename Rat_point_2>::iterator it = this->m_points.begin();
      if (this->m_points.size() == 1) {
        Rat_circle_2 circle(*it, m_bigints.back());
        Curve_2 conic(circle);
      } else {
        const Rat_point_2 & p1 = *it++;
        const Rat_point_2 & p2 = *it++;
        const Rat_point_2 & p3 = *it++;
        Curve_2 curve(p1, p2, p3);
        ++m_output_iterator = curve;
      }
    }
  };

  /*! Read the conic curves or arcs of conic curves  from the input file
   * \param filename the name of the input file
   * \param curves_out the iterator of the container of the read curves
   * \param bbox the counding box of the read curves
   */
  template<class OutputIterator>
  int read_data(const char * filename, OutputIterator curves_out,
                CGAL::Bbox_2 & bbox)
  {
    Conic_parser_visitor<OutputIterator> visitor(curves_out, bbox);
    if (!cb::benchmark_parse_file(filename, &visitor)) return -1;
    return 0;
  }

#if 0
  /*! */
  bool read_curve(std::ifstream & is, Curve_2 & cv)
  {
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
    } else if (type == 'c' || type == 'C') {
      // Construct a full circle. The line should have the format:
      //   c <x0> <y0> <R_sq>
      // where (x0, y0) is the center of the circle and R_sq is its squared
      // radius.
      CORE::BigInt    x0, y0, R_sq;
      
      str_line >> x0 >> y0 >> R_sq;
      Rat_point_2   p0(x0, y0);
      Rat_circle_2  circ(p0, R_sq);
      cv = Curve_2 (circ);
    } else if (type == 'e' || type == 'E') {
      CORE::BigInt    a, b, x0, y0;
      str_line >> a >> b >> x0 >> y0;

      Rational        r = b * b;
      Rational        s = a * a;
      Rational        u = -2 * r * x0;
      Rational        v = -2 * s * y0;
      Rational        w = r * x0 * x0 + s * y0 * y0 - r * s;
      cv = Curve_2 (r, s, 0, u, v, w);
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
#endif
};

#endif
