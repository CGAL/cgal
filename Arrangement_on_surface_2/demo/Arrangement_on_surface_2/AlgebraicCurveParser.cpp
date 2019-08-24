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

#include "AlgebraicCurveParser.h"

AlgebraicCurveParser::AlgebraicCurveParser(std::string &expression) : expression(expression) {};

bool AlgebraicCurveParser::validateExpression(const std::string &expression) {
  std::string expressionMutable = expression;
  for (auto iterator = expressionMutable.begin(); iterator != expressionMutable.end(); iterator++) {
    if (*iterator == 'x' || *iterator == 'y') continue;
    if (*iterator == '+' || *iterator == '-' || *iterator == '*' || *iterator == '^') continue;
    if (isdigit(*iterator) || isspace(*iterator)) continue;
    else return false;
  }
  return true;
}

Terms AlgebraicCurveParser::extractTerms() {
  Terms algebraicTerms;

  //Preprocess the expression: trim and remove spaces
  if (!this->validateExpression(expression)) return algebraicTerms;
  this->expression.erase(remove_if(expression.begin(), expression.end(), isspace), this->expression.end());


  using boost::spirit::ascii::space;
  typedef std::string::const_iterator iterator_type;
  typedef polynomial_parser<iterator_type> polynomial_parser;

  polynomial_parser pparser;
  std::string::const_iterator iter = expression.begin();
  std::string::const_iterator end = expression.end();

  //parsing goes on here
  bool r = phrase_parse(iter, end, pparser, space, algebraicTerms);

  if (r && iter == end)
  {
    std::cout<< "Parsing Successful";
    return algebraicTerms;
  }
  else
  {
    std::cout<<"Parsing Failed";
    return algebraicTerms;
  }
}
