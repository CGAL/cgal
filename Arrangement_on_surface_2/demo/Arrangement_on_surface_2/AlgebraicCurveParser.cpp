//
// Created by saurabh on 8/7/19.
//

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

  long termIndex = 0;

  //iterate through the string
  for (auto iterator = expression.begin(); iterator != expression.end(); iterator++) {
    auto next = iterator + 1;
    if (*next == '+' || *next == '-' || next == expression.end()) {
      //find the term between the +/- and the previous +/-
      std::string subExpression = std::string(expression.begin() + termIndex, next);
      AlgebraicCurveTerm term;
      bool result = parseTerm(subExpression.begin(), subExpression.end(), term);
      algebraicTerms.push_back(term);
      termIndex = next - expression.begin();
    }
  }
}

template<typename Iterator>
bool AlgebraicCurveParser::parseTerm(Iterator first, Iterator last, AlgebraicCurveTerm &term) {
  using boost::spirit::qi::long_long;
  using boost::spirit::qi::phrase_parse;
  using boost::phoenix::ref;
  long long xExponent = 0;
  long long yExponent = 0;
  long long coefficient = 1;
  phrase_parse(first, last,
          //Begin parser
               (
               long_long[ref(coefficient)] = boost::spirit::qi::_1 >>
                 -("x^" >> long_long[ref(
                 xExponent)] = boost::spirit::qi::_0)
                 >>
                 -("y^" >> long_long[ref(
                 yExponent)] = boost::spirit::qi::_0)
               |
                 "x^" >> long_long[ref(
                 xExponent)] = boost::spirit::qi::_0 >>
                 -("y^"
                 >> long_long[ref(
                                                                                                                   yExponent)] = boost::spirit::qi::_0)
               |
                  "y^" >> long_long[ref(
                  yExponent)] =
                  boost::spirit::qi::_0 >>
                  -("x^"
                  >> long_long[ref(
                  xExponent)] = boost::spirit::qi::_0)
               ),
               boost::spirit::ascii::space
  );
  return false;
}
