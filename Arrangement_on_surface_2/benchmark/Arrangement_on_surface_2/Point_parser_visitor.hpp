/*! \file
 * Parses points in a bff file
 */

#ifndef CURVE_PARSER_VISITOR_HPP
#define CURVE_PARSER_VISITOR_HPP

#include <CGAL/basic.h>
#include <CGAL/Bbox_2.h>
#include <string>

#include "CGAL/Benchmark/config.hpp"
#include "CGAL/Benchmark/Benchmark_visitor.hpp"
#include "CGAL/Benchmark/benchmark_format.hpp"
#include <CGAL/number_utils_classes.h>

#include "lexical_cast.hpp"
#include "Number_type_traits.hpp"

namespace cb = CGAL::benchmark;

/*! A portion of a visitor of the parser that reads points and constructs
 * their bounding box
 */
template <typename Rat_kernel, typename Rat_point, typename Rat_number_type>
class Point_parser_visitor : public cb::Benchmark_visitor {
protected:
  typedef std::list<Rat_point>                  Rat_point_list;
  typedef typename Rat_point_list::iterator     Rat_point_iter;

  Rat_point_list m_points;
  CGAL::Bbox_2 & m_bbox;

public:
  /*! Constructor */
  Point_parser_visitor(CGAL::Bbox_2 & bbox) : m_bbox(bbox) {}

  /*! Suppress output */
  virtual void token_not_handled(std::string s) {}

  virtual void accept_benchmark_name(std::string s)
  { Benchmark_visitor::accept_benchmark_name(s); }

  /*! A point maker */
  template <class FT>
  void make_point(FT & x, FT & y, Rat_point & p)
  {
    typename Rat_kernel::FT kx(x);
    typename Rat_kernel::FT ky(y);
    p = Rat_point(kx, ky);
    double dx = CGAL::to_double(x);
    double dy = CGAL::to_double(y);
    CGAL::Bbox_2 b(dx, dy, dx, dy);
    m_bbox = m_bbox + b;
  }

  /*! A point maker */
  template <class FT>
  void make_point(FT & x, FT & y, FT & w, Rat_point & p)
  {
    typename Rat_kernel::FT kx(x/w);
    typename Rat_kernel::FT ky(y/w);
    p = Rat_point(kx, ky);
  }

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
    typedef typename Number_type_traits<Rat_number_type>::FT            FT;
    FT x_ft = lexical_cast<FT>(x);
    FT y_ft = lexical_cast<FT>(y);
    Rat_point point;
    make_point(x_ft, y_ft, point);
    if (m_points.empty() || (point != m_points.back()))
      m_points.push_back(point);
    else
      std::cout << "Duplicate point: " << point << std::endl;
  }

  /*! Parse a generic Homogenuous point */
  virtual void accept_point_2( std::string x, std::string y, std::string w)
  {
    typedef typename Number_type_traits<Rat_number_type>::FT            FT;
    FT x_ft = lexical_cast<FT>(x);
    FT y_ft = lexical_cast<FT>(y);
    FT w_ft = lexical_cast<FT>(w);
    Rat_point point;
    make_point(x_ft, y_ft, point);
    if (m_points.empty() || (point != m_points.back()))
      m_points.push_back(point);
  }

  /*! Parse a rational Cartesian point */
  virtual void accept_point_2(std::string x_num, std::string x_denom,
                              std::string y_num, std::string y_denom)
  {
    typedef typename Number_type_traits<Rat_number_type>::RT            RT;
    typedef typename Number_type_traits<Rat_number_type>::FT            FT;
    typedef typename Number_type_traits<Rat_number_type>::Is_rational   Is_rational;
    RT x_num_rt = lexical_cast<RT>(x_num);
    RT x_denom_rt = lexical_cast<RT>(x_denom);
    RT y_num_rt = lexical_cast<RT>(y_num);
    RT y_denom_rt = lexical_cast<RT>(y_denom);
    FT x_ft, y_ft;
    make_ft(x_num_rt, x_denom_rt, x_ft, Is_rational());
    make_ft(y_num_rt, y_denom_rt, y_ft, Is_rational());
    Rat_point point;
    make_point(x_ft, y_ft, point);
    if (m_points.empty() || (point != m_points.back()))
      m_points.push_back(point);
  }
};

#endif
