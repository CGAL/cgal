// Copyright (c) 2012, 2020 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// SPDX-License-Identifier: GPL-3.0+
//
// Author(s): Saurabh Singh <ssingh@cs.iitr.ac.in>
//            Ahmed Essam <theartful.ae@gmail.com>

#include <CGAL/Polynomial_traits_d.h>
#include <CGAL/polynomial_utils.h>

#include <boost/spirit/include/phoenix.hpp>
#include <boost/spirit/include/qi.hpp>

#include "AlgebraicCurveParser.h"

namespace phx = boost::phoenix;
namespace qi = boost::spirit::qi;
namespace ascii = boost::spirit::ascii;

template <typename Polynomial_d, typename Iterator, typename Skipper>
struct PolynomialParser : qi::grammar<Iterator, Polynomial_d(), Skipper>
{
  using Self = PolynomialParser<Polynomial_d, Iterator, Skipper>;
  using Traits = CGAL::Polynomial_traits_d<Polynomial_d>;
  using Coefficient = typename Traits::Innermost_coefficient_type;
  using Innermost_leading_coefficient =
    typename Traits::Innermost_leading_coefficient;
  using Total_degree = typename Traits::Total_degree;

  PolynomialParser() : PolynomialParser::base_type(start)
  {
    using qi::_val;
    using qi::eps;

    for (int i = 0; i < Traits::d; i++)
        vars[i] = CGAL::shift(Polynomial_d(1), 1, i);

    // { expr = expr } or { expr }
    start = (expr >> '=' >> expr)[_val = qi::_1 - qi::_2] | expr[_val = qi::_1];
    // addition and subtraction
    expr = term[_val = qi::_1] >>
           *('+' >> term[_val += qi::_1] | '-' >> term[_val -= qi::_1]);
    // multiplication using *, and implied multiplication as in (x+y)(x+y)
    term = factor[_val = qi::_1] >>
           *(('*' >> factor[_val *= qi::_1]) | pow_expr[_val *= qi::_1]);
    // uniary - and + operators
    factor = qi::char_('-') >> pow_expr[_val = -qi::_1] |
             -qi::char_('+') >> pow_expr[_val = qi::_1];
    // power
    pow_expr =
      factor2[_val = qi::_1] >>
      (('^' >> factor2[_val = phx::bind(&Self::raise, this, _val, qi::_1)]) |
       eps);
    // ( expr )
    factor2 = (('(' >> expr >> ')') | factor3)[_val = qi::_1];
    // coefficients and variables
    factor3 =
      coeff[_val = phx::construct<Polynomial_d>(qi::_1)] | var[_val = qi::_1];
    coeff = qi::as_string[qi::lexeme[+qi::digit]]
                         [_val = phx::construct<Coefficient>(qi::_1)];
    if (Traits::d == 1)
      var = qi::char_('x')[_val = vars[0]];
    else
      var = qi::char_('x')[_val = vars[0]] | qi::char_('y')[_val = vars[1]];
  }

  Polynomial_d raise(const Polynomial_d& poly, const Polynomial_d& power)
  {
    if (total_degree(power) != 0)
    {
      this->error = true;
      return {};
    }

    return CGAL::ipower(
      poly, std::lround(CGAL::to_double(innermost_leading_coefficient(power))));
  }

  Innermost_leading_coefficient innermost_leading_coefficient;
  Total_degree total_degree;

  Polynomial_d vars[Traits::d];

  qi::rule<Iterator, Polynomial_d(), Skipper> start;
  qi::rule<Iterator, Polynomial_d(), Skipper> expr;
  qi::rule<Iterator, Polynomial_d(), Skipper> term;
  qi::rule<Iterator, Polynomial_d(), Skipper> pow_expr;
  qi::rule<Iterator, Polynomial_d(), Skipper> factor;
  qi::rule<Iterator, Polynomial_d(), Skipper> factor2;
  qi::rule<Iterator, Polynomial_d(), Skipper> factor3;
  qi::rule<Iterator, Polynomial_d(), Skipper> var;
  qi::rule<Iterator, Coefficient(), Skipper> coeff;

  bool error = false;
};

static bool hasValidChars2D(const std::string& expression)
{
  const char valid_chars[] = {'x', 'y', '+', '-', '*', '(', ')', '^', '='};
  return std::all_of(expression.begin(), expression.end(), [&](char c) {
    return std::isspace(c) || std::isdigit(c) ||
           std::find(std::begin(valid_chars), std::end(valid_chars), c) !=
             std::end(valid_chars);
  });
}

static bool hasValidChars1D(const std::string& expression)
{
  const char valid_chars[] = {'x', '+', '-', '*', '(', ')', '^'};
  return std::all_of(expression.begin(), expression.end(), [&](char c) {
    return std::isspace(c) || std::isdigit(c) ||
           std::find(std::begin(valid_chars), std::end(valid_chars), c) !=
             std::end(valid_chars);
  });
}

static inline bool hasValidChars(const std::string& expression, int dimension)
{
  if (dimension == 1)
    return hasValidChars1D(expression);
  else
    return hasValidChars2D(expression);
}

template <typename Polynomial_d>
boost::optional<Polynomial_d>
AlgebraicCurveParser<Polynomial_d>::operator()(const std::string& expression)
{
  using Traits = CGAL::Polynomial_traits_d<Polynomial_d>;
  using iterator_type = std::string::const_iterator;

  if (!hasValidChars(expression, Traits::d)) return {};

  PolynomialParser<Polynomial_d, iterator_type, ascii::space_type> pparser;
  std::string::const_iterator iter = expression.begin();
  std::string::const_iterator end = expression.end();

  // parsing goes on here
  Polynomial_d poly;
  bool r = qi::phrase_parse(iter, end, pparser, ascii::space, poly);

  if (r && iter == end && !pparser.error)
    return poly;
  else
    return {};
}

#ifdef CGAL_USE_CORE
// don't want to include ArrangementTypes.h
// makes compilation slower
// template class
// AlgebraicCurveParser<demo_types::Alg_seg_traits::Polynomial_2>;
// AlgebraicCurveParser<demo_types::Rational_traits::Polynomial_1>;
#include <CGAL/Algebraic_kernel_d_1.h>
#include <CGAL/Arr_algebraic_segment_traits_2.h>

template struct AlgebraicCurveParser<
  CGAL::Arr_algebraic_segment_traits_2<CORE::BigInt>::Polynomial_2>;

template struct AlgebraicCurveParser<
  typename CGAL::Algebraic_kernel_d_1<CORE::BigInt>::Polynomial_1>;
#endif
