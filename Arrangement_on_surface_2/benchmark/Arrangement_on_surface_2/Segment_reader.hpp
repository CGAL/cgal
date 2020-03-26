#ifndef SEGMENT_READER_HPP
#define SEGMENT_READER_HPP

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

/*! A reader of line segments */
template <class Traits>
class Segment_reader {
public:
  typedef typename Traits::Point_2              Point_2;
  typedef typename Traits::Curve_2              Curve_2;

  /*! A visitor of the parser that reads segements */
  template <class OutputIterator>
  class Segment_parser_visitor :
    public Point_parser_visitor<Kernel, Point_2, Number_type> {
  private:
    /*! The iterator of the output container */
    OutputIterator & m_output_iterator;

  public:
    /*! Constructor */
    Segment_parser_visitor(OutputIterator & oi, CGAL::Bbox_2 & bbox) :
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
      if (geom != "Lines")
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
      Curve_2 seg(this->m_points.front(), this->m_points.back());
      ++m_output_iterator = seg;
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
    Segment_parser_visitor<OutputIterator> visitor(curves_out, bbox);
    if (!cb::benchmark_parse_file(filename, &visitor)) return -1;
    return 0;
  }
};

#endif
