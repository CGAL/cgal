#ifndef SEGMENT_READER_H
#define SEGMENT_READER_H

#include <CGAL/basic.h>
#include <CGAL/Bbox_2.h>
#include <string>

#include "CGAL/benchmark_basic.h"
#include "CGAL/Benchmark_visitor.h"

#include "number_type.h"
#include "lexical_cast.h"
#include "Number_type_traits.h"
#include "Option_parser.h"

namespace cb = CGAL::benchmark;

template <class Traits>
class Segment_reader
{
public:
  typedef typename Traits::Point_2              Point_2;
  typedef typename Traits::Curve_2              Curve_2;

  /*! A visitor of the parser that reads segements */
  template <class OutputIterator>
  class Segment_parser_visitor : public cb::Benchmark_visitor {
  private:
    OutputIterator & m_output_iterator;
    Point_2 m_source;
    bool m_source_read;

    /*! Add a curve to the container the output iterator refers to
     * \param target the target point
     */
    void add_curve(Point_2 & target)
    {
      Curve_2 seg(m_source, target);
      ++m_output_iterator = seg;
    }
      
  public:
    /*! Constructor */
    Segment_parser_visitor(OutputIterator & oi) :
      m_output_iterator(oi),
      m_source_read(false)
    {}

    /*! Suppress output */
    virtual void token_not_handled(std::string s) {}

    virtual void accept_benchmark_name(std::string s)
    {
      Benchmark_visitor::accept_benchmark_name(s);
      std::cerr << "name '" << s << "', ";
    }

    /*! Accept only unbounded lines for Arrangements */ 
    virtual void accept_classification(std::string problem,
                                       std::string geom, 
                                       std::string clas,
                                       std::string family,
                                       std::string instance,
                                       std::string release)
    {
      if (problem != "Arrangement") error_handler( "classification error");
      if (geom != "Lines") error_handler( "classification error" );
      if (clas != "BoundedArcs") error_handler( "classification error" );
    }

    /*! A point maker */
    template <class FT>
    void make_point(FT & x, FT & y, Point_2 & p) { p = Point_2(x,y); }

    /*! A point maker */
    template <class FT>
    void make_point(FT & x, FT & y, FT & w, Point_2 & p)
    { p = Point_2(x/w, y/w); }

    /*! A point maker */
    template <class RT, class FT>
    void make_ft(RT & num, RT & denom, FT & n, CGAL::Tag_true)
    { n = FT(num, denom); }

    template <class RT, class FT>
    void make_ft(RT & num, RT & denom, FT & n, CGAL::Tag_false)
    { n = num / denom; }
    
    /*! Parse a generic Cartesian point */
    virtual void accept_point_2(std::string x, std::string y)
    {
      typedef Number_type_traits<Number_type>::FT       FT;
      FT x_ft = lexical_cast<FT>(x);
      FT y_ft = lexical_cast<FT>(y);
      if (m_source_read) {
        m_source_read = false;
        Point_2 p;
        make_point(x_ft, y_ft, p);
        add_curve(p);
        return;
      }
      make_point(x_ft, y_ft, m_source);
      m_source_read = true;
    }

    /*! Parse a generic Homogenuous point */
    virtual void accept_point_2( std::string x, std::string y, std::string w)
    {
      typedef Number_type_traits<Number_type>::FT       FT;
      FT x_ft = lexical_cast<FT>(x);
      FT y_ft = lexical_cast<FT>(y);
      FT w_ft = lexical_cast<FT>(w);
      if (m_source_read) {
        m_source_read = false;
        Point_2 p;
        make_point(x_ft, y_ft, p);
        add_curve(p);
        return;
      }
      make_point(x_ft, y_ft, w_ft, m_source);
      m_source_read = true;
    }

    /*! Parse a rational Cartesian point */
    virtual void accept_point_2(std::string x_num, std::string x_denom, 
                                std::string y_num, std::string y_denom)
    {
      typedef Number_type_traits<Number_type>::RT               RT;
      typedef Number_type_traits<Number_type>::FT               FT;
      typedef Number_type_traits<Number_type>::Is_rational      Is_rational;
      RT x_num_rt = lexical_cast<RT>(x_num);
      RT x_denom_rt = lexical_cast<RT>(x_denom);
      RT y_num_rt = lexical_cast<RT>(y_num);
      RT y_denom_rt = lexical_cast<RT>(y_denom);
      FT x_ft, y_ft;
      make_ft(x_num_rt, x_denom_rt, x_ft, Is_rational());
      make_ft(y_num_rt, y_denom_rt, y_ft, Is_rational());
      if (m_source_read) {
        m_source_read = false;
        Point_2 p;
        make_point(x_ft, y_ft, p);
        add_curve(p);
        return;
      }
      make_point(x_ft, y_ft, m_source);
      m_source_read = true;
    }
  };

  /*! Read the segments from the input file */
  template<class OutputIterator>
  int read_data(const char * filename, OutputIterator curves_out,
                CGAL::Bbox_2 & bbox)
  {
    Segment_parser_visitor<OutputIterator> visitor(curves_out);
    if (!cb::benchmark_parse_file(filename, &visitor)) {
      return -1;
    }
    return 0;
  }
};

#endif
