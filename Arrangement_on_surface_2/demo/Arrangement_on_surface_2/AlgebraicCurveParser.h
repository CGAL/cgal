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

#include <vector>
#include <string>
#include <complex>
#include <boost/config/warning_disable.hpp>
#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>
#include <boost/spirit/include/phoenix_object.hpp>
#include <boost/fusion/include/adapt_struct.hpp>
#include <boost/fusion/include/io.hpp>
#include <boost/bind.hpp>
#include <boost/spirit/include/phoenix.hpp>
#ifndef ARRANGEMENT_ON_SURFACE_2_DEMO_ALGEBRAICCURVEPARSERNEW_H
#define ARRANGEMENT_ON_SURFACE_2_DEMO_ALGEBRAICCURVEPARSERNEW_H
namespace qi = boost::spirit::qi;
namespace ascii = boost::spirit::ascii;

typedef std::vector<struct AlgebraicCurveTerm> Terms;
struct AlgebraicCurveTerm {
    boost::optional<long long> xExponent, yExponent, coefficient;
};

template <typename Iterator>
    struct polynomial_parser : qi::grammar<Iterator,
        std::vector<AlgebraicCurveTerm>(),
        ascii::space_type>
    {
        polynomial_parser() : polynomial_parser::base_type(start)
        {
            namespace phx = boost::phoenix;
            using qi::int_;
            using qi::lit;
            using qi::double_;
            using qi::lexeme;
            using ascii::char_;
            using qi::_val;
            using qi::eps;

            exponent %= '^' >> int_;

            xExponent %= 'x' >> ( exponent | eps[_val = 1] );

            yExponent %= 'y' >> ( exponent | eps[_val = 1] );

            poly_term = eps[_val = AlgebraicCurveTerm()]
                  >>
                    -int_[phx::bind(&AlgebraicCurveTerm::coefficient, _val) = qi::_1]
                  >>
                    -xExponent[phx::bind(&AlgebraicCurveTerm::xExponent, _val) = qi::_1]
                  >>
                    -yExponent[phx::bind(&AlgebraicCurveTerm::yExponent, _val) = qi::_1]
            ;

            negative_poly_term = eps[_val = AlgebraicCurveTerm()] >>
                (
                    int_[phx::bind(&AlgebraicCurveTerm::coefficient, _val) = -1 * qi::_1]
                    | eps[phx::bind(&AlgebraicCurveTerm::coefficient, _val) = -1]
                ) >>
                    -xExponent[phx::bind(&AlgebraicCurveTerm::xExponent, _val) = qi::_1]
                  >>
                    -yExponent[phx::bind(&AlgebraicCurveTerm::yExponent, _val) = qi::_1]
            ;

            start = eps[_val = std::vector<AlgebraicCurveTerm>()]
                >> poly_term[phx::push_back(_val, qi::_1)]
                >> *(
                    ('+' >> poly_term[phx::push_back(_val, qi::_1)])
                    |
                    ('-' >> negative_poly_term[phx::push_back(_val, qi::_1)])
                    )
            ;
        }

        qi::rule<Iterator, int(), ascii::space_type> exponent;
        qi::rule<Iterator, int(), ascii::space_type> xExponent;
        qi::rule<Iterator, int(), ascii::space_type> yExponent;
        qi::rule<Iterator, AlgebraicCurveTerm(), ascii::space_type> poly_term;
        qi::rule<Iterator, AlgebraicCurveTerm(), ascii::space_type> negative_poly_term;
        qi::rule<Iterator, std::vector<AlgebraicCurveTerm>(), ascii::space_type> start;
    };

class AlgebraicCurveParser {
public:
    explicit AlgebraicCurveParser(std::string& expression);
    bool validateExpression(const std::string &expression);
    Terms extractTerms();
    std::string expression;
};


#endif //ARRANGEMENT_ON_SURFACE_2_DEMO_ALGEBRAICCURVEPARSERNEW_H
