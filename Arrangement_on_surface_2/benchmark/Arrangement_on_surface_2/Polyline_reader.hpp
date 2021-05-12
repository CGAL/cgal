#ifndef POLYLINE_READER_HPP
#define POLYLINE_READER_HPP

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

/*! Convert Number_type to int. Used to calculate the bounding box */
template <class T_NT>
class Toint {
public:
  int operator()(T_NT & n) { return static_cast<int>(n); }
};

/*! Convert leda_rational to int. Used to calculate the bounding box */
#if BENCH_NT == LEDA_RAT_NT
template <>
class Toint<leda_rational> {
public:
  int operator()(leda_rational & n)
  {
    if (n.denominator() == n.numerator()) return 1;
    leda_integer lint = n.denominator() / n.denominator();
    return lint.to_long();
  }
};
#endif

/*! A reader of polyline curves */
template <class Traits>
class Polyline_reader {
public:
  typedef typename Traits::Point_2              Point_2;
  typedef typename Traits::Curve_2              Curve_2;

  /*! A parser of bff files that contains polylines */
  template <typename OutputIterator>
  class Polyline_parser_visitor :
    public Point_parser_visitor<Kernel, Point_2, Number_type> {
  private:
    /*! The iterator of the output container */
    OutputIterator & m_output_iterator;

  public:
    /*! Constructor */
    Polyline_parser_visitor(OutputIterator & oi, CGAL::Bbox_2 & bbox) :
      Point_parser_visitor<Kernel, Point_2, Number_type>(bbox),
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
      if (geom != "Polylines")
        this->error_handler( "classification error" );
      if (clas != "BoundedArcs")
        this->error_handler( "classification error" );
    }

    /*! Process start polyline */
    virtual void begin_polyline_2() { this->m_points.clear(); }

    /*! Process end polyline
     * Add a curve to the container the output iterator refers to
     */
    virtual void end_polyline_2()
    {
      Curve_2 polyline(this->m_points.begin(), this->m_points.end());
      ++m_output_iterator = polyline;
    }
  };

  /*! Read the segments from the input file
   * \param filename the name of the input file
   * \param curves_out the iterator of the container of the read curves
   * \param bbox the counding box of the read curves
   */
  template<class OutputIterator>
  int read_data(const char * filename, OutputIterator curves_out,
                CGAL::Bbox_2 & bbox)
  {
    Polyline_parser_visitor<OutputIterator> visitor(curves_out, bbox);
    if (!cb::benchmark_parse_file(filename, &visitor)) return -1;
    return 0;
  }
};

#endif
