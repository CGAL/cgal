//
// Created by saurabh on 8/7/19.
//
#include <vector>
#include <string>
#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>
#ifndef ARRANGEMENT_ON_SURFACE_2_DEMO_ALGEBRAICCURVEPARSERNEW_H
#define ARRANGEMENT_ON_SURFACE_2_DEMO_ALGEBRAICCURVEPARSERNEW_H

typedef std::vector<struct AlgebraicCurveTerm> Terms;
struct AlgebraicCurveTerm {
    long long xExponent, yExponent, coefficient;
    bool isPositive;
};
class AlgebraicCurveParser {
public:
    explicit AlgebraicCurveParser(std::string& expression);
    bool validateExpression(const std::string &expression);
    Terms extractTerms();
    std::string expression;
private:
    template <typename Iterator>
    bool parseTerm(Iterator first, Iterator last, AlgebraicCurveTerm& term);
    bool extractSign(std::string subExpression);
    bool signPresent(std::string subExpression);


};


#endif //ARRANGEMENT_ON_SURFACE_2_DEMO_ALGEBRAICCURVEPARSERNEW_H
