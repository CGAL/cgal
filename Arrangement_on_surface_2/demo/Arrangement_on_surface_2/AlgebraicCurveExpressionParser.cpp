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
// Author(s): Tianyu Zhou <zhoutianyu16tue@gmail.com>

#include "AlgebraicCurveExpressionParser.h"

AlgebraicCurveExpressionParser::
AlgebraicCurveExpressionParser(std::string& poly_expr):
  poly_expr(poly_expr)
{}

//! preprocessing of the expression
/*!
  \param str string argument of the polynomial expression
*/
void AlgebraicCurveExpressionParser::pre_hanlde_poly_expr( std::string& str )
{
  std::vector<int> neg_sign_index;

  for (size_t i = 0; i < str.size(); ++i)
  {
    if ( str[i] == '+' )
    {
      str.replace(i, 1, ",");
    }

    if ( str[i] == '-')
    {
      neg_sign_index.push_back(i);
    }

    if ( str[i] == '*')
    {
      str.replace(i, 1, " ");
    }
  }

  for (int i = int(neg_sign_index.size()-1); i >=0 ; i--)
  {
    str.insert(neg_sign_index[i], "," );
  }
}

//! extracting the coefficients from the expression
/*!
  \param term structure to store coefficients and the exponents
  \param poly_expr polynomial expression itselfs
  \return integer values
*/
int AlgebraicCurveExpressionParser::
extract_poly_coefficient(std::string& poly_expr, struct term& term)
{
  int ret_val = 0;

  if ( poly_expr[0] == '-' && (poly_expr[1]=='x' || poly_expr[1]=='y') )
  {
    term.coefficient = -1;
    ret_val = 1;
  }
  else if (poly_expr[0]=='x' || poly_expr[0]=='y')
  {
    term.coefficient = 1;
    ret_val = 0;
  }
  else{
    term.coefficient = std::stoi(poly_expr);
    std::string tmp = std::to_string(term.coefficient);
    ret_val = int(tmp.size());
  }
  return ret_val;
}

//! extracting the exponents from the expression
/*!
  \param term structure to store coefficients and the exponents
  \param poly_expr polynomial expression itselfs
*/
void AlgebraicCurveExpressionParser::
extract_poly_exponents(std::string& sub_poly_expr, struct term& term)
{
  term.x_exponent = 0;
  term.y_exponent = 0;

  std::size_t found_x = sub_poly_expr.find('x');
  if (found_x != std::string::npos)
  {
    if (sub_poly_expr[found_x+1] == '^')
    {
      term.x_exponent = std::stoi(sub_poly_expr.substr(found_x+2));
    }
    else
    {
      term.x_exponent = 1;
    }
  }

  std::size_t found_y = sub_poly_expr.find('y');
  if (found_y != std::string::npos)
  {
    if (sub_poly_expr[found_y+1] == '^')
    {
      term.y_exponent = std::stoi(sub_poly_expr.substr(found_y+2));
    }
    else
    {
      term.y_exponent = 1;
    }
  }
}

//! printing out the expression correctly
/*!
  \param term structure that contains coefficients and exponents
*/
void AlgebraicCurveExpressionParser::print_poly_term( const term& term )
{
    std::cout << "term\n";
    std::cout << "  coefficient: " << term.coefficient << "\n"
              << "  x_exponent: " << term.x_exponent << "\n"
              << "  y_exponent: " << term.y_exponent << "\n";
}

//! function to handle correct ordering of extraction
/*!
  \param term structure that contains coefficients and the exponents
  \param poly_expr polynomial expression itselfs
*/
void AlgebraicCurveExpressionParser::
extract_poly_components(std::string& poly_expr, struct term& term)
{
  // Extract the coefficient
  int starting_index = extract_poly_coefficient(poly_expr, term);

  // Extract the x y exponents respectively
  std::string sub_poly_expr = poly_expr.substr(starting_index);
  extract_poly_exponents(sub_poly_expr, term);
}

//! grouping terms together
/*!
  \param poly_expr vector of structure term that holds the grouping
*/
bool AlgebraicCurveExpressionParser::
extract_poly_terms(std::vector<struct term>& poly_terms)
{
  // Remove all white spaces
  this->poly_expr.erase(remove_if(this->poly_expr.begin(),
                                  this->poly_expr.end(), isspace),
                        this->poly_expr.end());
  if (this->poly_expr.size() == 0)
  {
    return false;
  }

  boost::char_separator<char> equal_sign_sep("=");
  boost::tokenizer< boost::char_separator<char> >
    equation_sides(this->poly_expr, equal_sign_sep);

  std::vector<std::vector<struct term>> sides;

  BOOST_FOREACH (const std::string& str, equation_sides) {

    std::string equation_side = str;
    pre_hanlde_poly_expr(equation_side);

    boost::char_separator<char> sep(",");
    boost::tokenizer< boost::char_separator<char> > tokens(equation_side, sep);

    std::vector<struct term> terms_temp;

    BOOST_FOREACH (const std::string& str, tokens) {
      std::string output = str;

      struct term term;

      extract_poly_components(output, term);
      terms_temp.push_back(term);
    }

    sides.push_back(terms_temp);
  }

  if (sides.size() > 2)
  {
    return false;
  }

  std::vector<struct term> left_side = sides[0];
  std::vector<struct term> right_side;
  if (sides.size() == 2)
  {
    right_side = sides[1];
  }

  for (size_t i = 0; i < left_side.size(); ++i)
  {
    // print_poly_term(left_side[i]);
    poly_terms.push_back(left_side[i]);
  }

  for (size_t i = 0; i < right_side.size(); ++i)
  {
    right_side[i].coefficient = -right_side[i].coefficient;
    // print_poly_term(right_side[i]);
    poly_terms.push_back(right_side[i]);
  }

  return true;
}
