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

template <>
CORE::BigRat lexical_cast<CORE::BigRat>(std::string & str)
{
  std::stringstream ss(str);
  CORE::BigRat result;
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
    public Point_parser_visitor<Rat_kernel, Rat_point_2, Rational> {
  private:
    typedef Point_parser_visitor<Rat_kernel, Rat_point_2, Rational> Base;
    
    /*! The iterator of the output container */
    OutputIterator & m_output_iterator;

    /*! Big ints used to store coefficients */
    std::list<CORE::BigInt> m_bigints;

    /*! Is a conivs arc currently being processed? */
    bool m_processing_arc;

    /*! A place holder to store the undelying conic of a conic arc */
    Curve_2 m_conic;

    /*! Last orientation */
    CGAL::Orientation m_orient;

  public:
    /*! Constructor */
    Conic_parser_visitor(OutputIterator & oi, CGAL::Bbox_2 & bbox) :
      Point_parser_visitor<Rat_kernel, Rat_point_2, Rational>(bbox),
      m_output_iterator(oi),
      m_processing_arc(false),
      m_orient(CGAL::COUNTERCLOCKWISE)
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

    virtual void accept_integer(std::string s)
    {
      CORE::BigInt integer = lexical_cast<CORE::BigInt>(s);
      m_bigints.push_back(integer);
    }
    
    virtual void accept_orientation( std::string s)
    {
      m_orient =
        (s == "COUNTERCLOCKWISE") ? CGAL::COUNTERCLOCKWISE : CGAL::CLOCKWISE;
    }

    /*! Start a full curve */
    void begin_full()
    {
      this->m_points.clear();
      m_bigints.clear();
    }

    /*! Start an arc */
    void begin_arc()
    {
      begin_full();
      m_processing_arc = true;
    }
    
    /*! Start a line segment */
    virtual void begin_line_segment_2() { begin_arc(); }

    /*! End a line segment */
    virtual void end_line_segment_2()
    {
      Rat_segment_2 seg (this->m_points.front(), this->m_points.back());
      Curve_2 conic(seg);
      ++m_output_iterator = conic;
    }

    /*! Start a circle */
    virtual void begin_circle_2() { begin_full(); }

    /*! End a circle */
    virtual void end_circle_2()
    {
      if (this->m_points.size() == 1) {
        Rat_point_2 & point = this->m_points.back();
        Rational radius_square = m_bigints.back();
        Rat_circle_2 circle(point, radius_square);
        m_conic = circle;
        if (!m_processing_arc) ++m_output_iterator = m_conic;
        this->m_points.clear();

        double radius = sqrt(CGAL::to_double(radius_square));
        double x = CGAL::to_double(point.x());
        double y = CGAL::to_double(point.y());
        CGAL::Bbox_2 b(x - radius, y - radius, x + radius, y + radius);
        this->m_bbox = this->m_bbox + b;
        return;
      }

      typename Base::Rat_point_iter it = this->m_points.begin();
      const Rat_point_2 & p1 = *it++;
      const Rat_point_2 & p2 = *it++;
      const Rat_point_2 & p3 = *it++;
      Rat_circle_2 circle(p1, p2, p3);
      m_conic = circle;
      if (!m_processing_arc) ++m_output_iterator = m_conic;
      this->m_points.clear();
    }

    /*! End an arc */
    void end_arc()
    {
      m_processing_arc = false;
      CGAL_assertion(this->m_points.size() == 2);
      Rat_point_2 & ps_rat = this->m_points.front();
      Rat_point_2 & pt_rat = this->m_points.back();
      Point_2 ps(ps_rat.x(), ps_rat.y());
      Point_2 pt(pt_rat.x(), pt_rat.y());
      Curve_2 conic_arc(m_conic.r(), m_conic.s(), m_conic.t(),
                        m_conic.u(), m_conic.v(), m_conic.w(),
                        m_orient, ps ,pt);
      ++m_output_iterator = conic_arc;
    }
    
    /*! Start a circle arc */
    virtual void begin_circle_arc_2() { begin_arc(); }

    /*! End a circle arc */
    virtual void end_circle_arc_2()
    {
      if (this->m_points.size() == 3) {
        typename Base::Rat_point_iter it = this->m_points.begin();
        const Rat_point_2 & p1 = *it++;
        const Rat_point_2 & p2 = *it++;
        const Rat_point_2 & p3 = *it++;
        Curve_2 curve(p1, p2, p3);
        ++m_output_iterator = curve;
        return;
      }
      end_arc();
    }

    /*! Start an iso-ellipse */
    virtual void begin_iso_ellipse_2() { begin_full(); }

    /*! End an iso-ellipse */
    virtual void end_iso_ellipse_2()
    {
      std::list<CORE::BigInt>::iterator ii = this->m_bigints.end();
      Rational r = *(--ii);
      Rational s = *(--ii);

      Rational x = this->m_points.back().x();
      Rational y = this->m_points.back().y();

      Rational u = -2 * r * x;
      Rational v = -2 * s * y;
      Rational w = r * x * x + s * y * y - r * s;

      m_conic = Curve_2(r, s, 0, u, v, w);

      double rx = sqrt(CGAL::to_double(r));
      double ry = sqrt(CGAL::to_double(s));
      double xd = CGAL::to_double(x);
      double yd = CGAL::to_double(y);
      CGAL::Bbox_2 b(xd - rx, yd - rx, xd + ry, yd + ry);
      this->m_bbox = this->m_bbox + b;

      if (!m_processing_arc) ++m_output_iterator = m_conic;
      this->m_points.clear();
    }

    /*! Start an ellipse arc */
    virtual void begin_iso_ellipse_arc_2() { begin_arc(); }

    /*! End an ellipse arc */
    virtual void end_iso_ellipse_arc_2() { end_arc(); }
    
    /*! Accept a conic curve */
    virtual void accept_conic_2(std::string r_str, std::string s_str,
                                std::string t_str, std::string u_str,
                                std::string v_str, std::string w_str)
    {
      CORE::BigInt r = lexical_cast<CORE::BigInt>(r_str);
      CORE::BigInt s = lexical_cast<CORE::BigInt>(s_str);
      CORE::BigInt t = lexical_cast<CORE::BigInt>(t_str);
      CORE::BigInt u = lexical_cast<CORE::BigInt>(u_str);
      CORE::BigInt v = lexical_cast<CORE::BigInt>(v_str);
      CORE::BigInt w = lexical_cast<CORE::BigInt>(w_str);
      m_conic = Curve_2(r, s, t, u, v, w);
      if (!m_processing_arc) ++m_output_iterator = m_conic;
    }

    /*! Start a conic arc */
    virtual void begin_conic_arc_2() { begin_arc(); }

    /*! End a conic arc */
    virtual void end_conic_arc_2() { end_arc(); }
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
  }
#endif
};

#endif
