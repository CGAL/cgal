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

#include "AlgebraicCurveParserOld.h"

AlgebraicCurveParserOld::AlgebraicCurveParserOld(std::string& expression) :
  expression(expression)
{};

TermsArray AlgebraicCurveParserOld::extractTerms() {
  TermsArray algebraicTerms;

  //Preprocess the expression: trim and remove spaces
  if (!this->validateExpression(expression)) return algebraicTerms;
  this->expression.erase(remove_if(expression.begin(), expression.end(),
                                   isspace), this->expression.end());

  long termIndex = 0;

  //iterate through the string
  for (auto iterator = expression.begin(); iterator != expression.end();
       iterator++)
  {
    auto next = iterator + 1;
    if (*next == '+' || *next == '-' || next == expression.end()) {
      //find the term between the +/- and the previous +/-
      std::string subExpression = std::string(expression.begin() + termIndex, next);
      AlgebraicTerm term = extractCoefficientAndExponent(subExpression);
      algebraicTerms.push_back(term);
      termIndex = next - expression.begin();
    }
  }
  return algebraicTerms;
}

AlgebraicTerm AlgebraicCurveParserOld::
extractCoefficientAndExponent(std::string& subExpression) {
  AlgebraicTerm algebraicTerm{};

  algebraicTerm.positiveSign = extractSign(subExpression);
  algebraicTerm.coefficient = extractCoefficient(subExpression);
  algebraicTerm.xExponent = extractXExponent(subExpression);
  algebraicTerm.yExponent = extractYExponent(subExpression);
  return algebraicTerm;
}

int AlgebraicCurveParserOld::extractXExponent(std::string& subExpression) {
  bool extractingExponent = false;
  int exponentValue = 0;
  std::string exponentString;

  for (auto iterator = subExpression.begin(); iterator != subExpression.end(); iterator++) {
    if (extractingExponent) {
      if (*iterator == '^') {
        continue;
      }
      else if (*iterator == 'y' || *iterator == '*') {
        extractingExponent = false;
        std::string::size_type sz;
        if (!exponentString.empty())
          exponentValue += std::stoi(exponentString, &sz);
        continue;
      }
      else {
        char charAtPosition = *iterator;
        exponentString += charAtPosition;
        if (iterator == subExpression.end() - 1) {
          extractingExponent = false;
          std::string::size_type sz;
          if (!exponentString.empty())
            exponentValue += std::stoi(exponentString, &sz);
          continue;
        }
      }
    }
    else {
      if (*iterator == 'x') {
        if (iterator == subExpression.end() - 1) {
          exponentValue += 1;
          continue;
        } else if (*(iterator + 1) != '^') {
          exponentValue += 1;
        } else if (*(iterator + 1) == '^') {
          extractingExponent = true;
          exponentString = "";
          continue;
        }
      }
    }
  }

  return exponentValue;
}

int AlgebraicCurveParserOld::extractYExponent(std::string& subExpression) {
  bool extractingExponent = false;
  int exponentValue = 0;
  std::string exponentString;

  for (auto iterator = subExpression.begin(); iterator != subExpression.end(); iterator++) {
    if (extractingExponent) {
      if (*iterator == '^') {
        continue;
      }
      else if (*iterator == 'y' || *iterator == '*') {
        extractingExponent = false;
        std::string::size_type sz;
        if (!exponentString.empty()) exponentValue += std::stoi(exponentString, &sz);
        continue;
      }
      else {
        char charAtPosition = *iterator;
        exponentString += charAtPosition;
        if (iterator == subExpression.end() - 1) {
          extractingExponent = false;
          std::string::size_type sz;
          if (!exponentString.empty()) exponentValue += std::stoi(exponentString, &sz);
          continue;
        }
      }
    }
    else {
      if (*iterator == 'y') {
        if (iterator == subExpression.end() - 1) {
          exponentValue += 1;
          continue;
        } else if (*(iterator + 1) != '^') {
          exponentValue += 1;
        } else if (*(iterator + 1) == '^') {
          extractingExponent = true;
          exponentString = "";
          continue;
        }
      }
    }
  }

  return exponentValue;
}

int AlgebraicCurveParserOld::extractCoefficient(std::string& subExpression) {
  bool signPresent = this->signPresent(subExpression);
  std::string coefficientString;
  int coefficientValue;

  std::string::iterator iteratorStart;
  if (signPresent) iteratorStart = subExpression.begin() + 1;
  else iteratorStart = subExpression.begin();
  if (*iteratorStart == 'x' || *iteratorStart == 'y') {
    return 1;
  }

  for (auto iterator = iteratorStart; iterator != subExpression.end(); iterator++) {
    //determine coefficient
    if (*iterator == 'x' || *iterator == 'y' || *iterator == '*') {
      break;
    }
    char value = *iterator;
    coefficientString += value;
  }

  std::string::size_type sz;
  coefficientValue = std::stoi(coefficientString, &sz);

  return coefficientValue;
}

bool AlgebraicCurveParserOld::extractSign(std::string& subExpression) {
  bool positiveSign;
  switch (*subExpression.begin()) {
   case '+':
    positiveSign = true;
    break;
   case '-':
    positiveSign = false;
    break;
   default:
    positiveSign = true;
  }
  return positiveSign;
}

bool AlgebraicCurveParserOld::signPresent(std::string& subExpression) {
  bool signPresent;
  switch (*subExpression.begin()) {
   case '+':
    signPresent = true;
    break;
   case '-':
    signPresent = true;
    break;
   default:
    signPresent = false;
  }
  return signPresent;
}

bool AlgebraicCurveParserOld::validateExpression(const std::string& expression) {
  //check for unwanted variables:
  std::string expressionMutable = expression;
  for (auto iterator = expressionMutable.begin(); iterator!=expressionMutable.end(); iterator++)
  {
    if (*iterator == 'x' || *iterator == 'y') continue;
    if (*iterator == '+' || *iterator =='-' || *iterator == '*' || *iterator == '^' ) continue;
    if (isdigit(*iterator) || isspace(*iterator)) continue;
    else return false;
  }
  return true;
}
