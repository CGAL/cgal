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
// Author(s)     : Saurabh Singh <ssingh@cs.iitr.ac.in>
//
#include <string>
#include <vector>
#include <algorithm>

#ifndef ARRANGEMENT_ON_SURFACE_2_DEMO_ALGEBRAICCURVEPARSER_H
#define ARRANGEMENT_ON_SURFACE_2_DEMO_ALGEBRAICCURVEPARSER_H

typedef std::vector<struct AlgebraicTerm> TermsArray;

struct AlgebraicTerm {
    int xExponent;
    int yExponent;
    int coefficient;
    bool positiveSign;
};

class AlgebraicCurveParserOld {
public:
    AlgebraicCurveParserOld(std::string &expression);

    bool validateExpression(const std::string &expression);

    TermsArray extractTerms();

    std::string expression;

private:

    AlgebraicTerm extractCoefficientAndExponent(std::string &subExpression);

    int extractXExponent(std::string &subExpression);

    int extractYExponent(std::string &subExpression);

    int extractCoefficient(std::string &subExpression);

    bool extractSign(std::string &subExpression);

    bool signPresent(std::string &subExpression);
};


#endif //ARRANGEMENT_ON_SURFACE_2_DEMO_ALGEBRAICCURVEPARSER_H
