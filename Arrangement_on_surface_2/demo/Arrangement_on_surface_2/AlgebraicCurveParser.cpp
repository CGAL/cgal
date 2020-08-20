//
// Copyright (c) 2012  Tel-Aviv University (Israel).
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

#include <CGAL/Arr_algebraic_segment_traits_2.h>
#include <CGAL/Polynomial_traits_d.h>

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/phoenix.hpp>

#include "AlgebraicCurveParser.h"

namespace phx = boost::phoenix;
namespace qi = boost::spirit::qi;
namespace ascii = boost::spirit::ascii;

template <typename Polynomial_2, typename Iterator>
struct PolynomialParser :
    qi::grammar<Iterator, Polynomial_2(), ascii::space_type>
{
  using Traits = CGAL::Polynomial_traits_d<Polynomial_2>;
  using Coefficient = typename Traits::Innermost_coefficient_type;

  PolynomialParser() : PolynomialParser::base_type(start)
  {
    using ascii::space;
    using qi::_val;
    using qi::eps;

    auto raise = [&](auto&& poly, auto&& power) {
      if (power.degree() != 0)
      {
        error = true;
        return Polynomial_2{};
      }
      return CGAL::ipower(poly, power.begin()->begin()->intValue());
    };

    // expr = expr or expr
    start = (expr >> '=' >> expr)[_val = qi::_1 - qi::_2] | expr[_val = qi::_1];
    // addition and subtraction
    expr = term[_val = qi::_1] >>
           *('+' >> term[_val += qi::_1] | '-' >> term[_val -= qi::_1]);
    // multiplication using *, and implied multiplication as in (x+y)(x+y)
    term = pow_expr[_val = qi::_1] >>
           *((qi::char_('*') >> pow_expr[_val *= qi::_1]) |
             pow_expr2[_val *= qi::_1]);
    // power
    pow_expr = factor[_val = qi::_1] >>
               (('^' >> factor2[_val = phx::bind(raise, _val, qi::_1)]) | eps);
    // to avoid ambiguity in a case like x^2+y
    // which can be seen as implied multiplication (x^2)*(+y)
    pow_expr2 = factor2[_val = qi::_1] >>
                (('^' >> factor2[_val = phx::bind(raise, _val, qi::_1)]) | eps);
    // uniary - and + operators
    factor = qi::char_('-') >> factor2[_val = -qi::_1] |
             -qi::char_('+') >> factor2[_val = qi::_1];
    // ( expr )
    factor2 = (('(' >> expr >> ')') | factor3)[_val = qi::_1];
    // coefficients and variables
    factor3 =
      coeff[_val = phx::construct<Polynomial_2>(qi::_1)] | var[_val = qi::_1];
    coeff = qi::as_string[qi::lexeme[+qi::digit]]
                         [_val = phx::construct<Coefficient>(qi::_1)];
    var = qi::char_('x')[_val = x] | qi::char_('y')[_val = y];
  }

  Polynomial_2 x = CGAL::shift(Polynomial_2(1), 1, 0);
  Polynomial_2 y = CGAL::shift(Polynomial_2(1), 1, 1);
  qi::rule<Iterator, Polynomial_2(), ascii::space_type> start;
  qi::rule<Iterator, Polynomial_2(), ascii::space_type> expr;
  qi::rule<Iterator, Polynomial_2(), ascii::space_type> term;
  qi::rule<Iterator, Polynomial_2(), ascii::space_type> pow_expr;
  qi::rule<Iterator, Polynomial_2(), ascii::space_type> pow_expr2;
  qi::rule<Iterator, Polynomial_2(), ascii::space_type> factor;
  qi::rule<Iterator, Polynomial_2(), ascii::space_type> factor2;
  qi::rule<Iterator, Polynomial_2(), ascii::space_type> factor3;
  qi::rule<Iterator, Polynomial_2(), ascii::space_type> var;
  qi::rule<Iterator, Coefficient(), ascii::space_type> coeff;
  bool error = false;
};

static inline bool hasValidChars(const std::string& expression)
{
  const char valid_chars[] = {'x', 'y', '+', '-', '*', '(', ')', '^', '='};
  return std::all_of(expression.begin(), expression.end(), [&](char c) {
    return std::isspace(c) || std::isdigit(c) ||
           std::find(std::begin(valid_chars), std::end(valid_chars), c) !=
             std::end(valid_chars);
  });
}

template <typename Polynomial_2>
boost::optional<Polynomial_2>
AlgebraicCurveParser<Polynomial_2>::operator()(const std::string& expression)
{
  using ascii::space;
  using iterator_type = std::string::const_iterator;

  if (!hasValidChars(expression)) return {};

  PolynomialParser<Polynomial_2, iterator_type> pparser;
  std::string::const_iterator iter = expression.begin();
  std::string::const_iterator end = expression.end();

  // parsing goes on here
  Polynomial_2 poly;
  bool r = qi::phrase_parse(iter, end, pparser, space, poly);

  if (r && iter == end && !pparser.error)
    return poly;
  else
    return {};
}

// don't want to include ArrangementTypes.h
// makes compilation slower
// template class
// AlgebraicCurveParser<demo_types::Alg_seg_traits::Polynomial_2>;
template struct AlgebraicCurveParser<
  CGAL::Arr_algebraic_segment_traits_2<CORE::BigInt>::Polynomial_2>;
