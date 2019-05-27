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
// Author(s)     : Tianyu Zhou <zhoutianyu16tue@gmail.com>


#ifndef ALGEBRAIC_CURVE_EXPRESSION_PARSER
#define ALGEBRAIC_CURVE_EXPRESSION_PARSER

#include <iostream>
#include <string>
#include <vector>
#include <boost/foreach.hpp>
#include <boost/tokenizer.hpp>

struct term
{
    int coefficient;
    int x_exponent;
    int y_exponent;
};

class AlgebraicCurveExpressionParser
{
public:
    AlgebraicCurveExpressionParser(std::string& poly_expr);

    bool extract_poly_terms(std::vector<struct term>& poly_terms);                  //!< extracting the terms seperately 
    
private:
    void pre_hanlde_poly_expr( std::string& str );                                  //!< preprocessing of the polynomial
    int extract_poly_coefficient(std::string& poly_expr, struct term& term);        //!< extracting the coefficients from the polynomila
    void extract_poly_exponents(std::string& sub_poly_expr, struct term& term);     //!< extracting exponents from the polynomial
    void print_poly_term( const term& term );                                       //!< debugger to see if the extraction works
    void extract_poly_components(std::string& poly_expr, struct term& term);

    std::string poly_expr;                                                          /*!< polynomial expression of string */
};

#endif // ALGEBRAIC_CURVE_EXPRESSION_PARSER
