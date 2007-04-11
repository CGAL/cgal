/**************************************************************************
// Copyright (c) 2004-2007  Max-Planck-Institut Saarbruecken (Germany)
// All rights reserved.
//
// This file is part of BenchmarkParser; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with BenchmarkParser.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
//
// Author(s) : Andreas Meyer
**************************************************************************/

#ifndef CGAL_BENCHMARK_POLYNOMIAL_VISITOR_HPP
#define CGAL_BENCHMARK_POLYNOMIAL_VISITOR_HPP

#include <string>
#include <vector>

#include <CGAL/Object.h>
#include <CGAL/Benchmark/benchmark_visitor.hpp>

CGAL_BENCHMARK_BEGIN_NAMESPACE

class Polynomial_benchmark_visitor
  : public Benchmark_visitor 
{
    typedef std::vector< Object >       Object_container;
    typedef Object_container::iterator  iterator;
    
    // feel free to extend this list.
    //enum Coefficient_numbertype {
    //  Not_supported = 0, 
    //  Integer, 
    //  Sqrt_ext_int_int 
    //};
    
    bool                        inside_monom;
    
    Object                      coefficient_value;
    std::vector< unsigned int > exponent_vector;
    
    unsigned int                number_of_variables;
    //Coefficient_numbertype      coefficient_numbertype;
    
    
    Object_container            object_container;
    
public:
    Polynomial_visitor();

    //virtual void token_not_handled( std::string s) {}

    // accept everything
    virtual void accept_classification( std::string problem,
                                        std::string geom,
                                        std::string clas,
                                        std::string family,
                                        std::string instance,
                                        std::string release) {}


    virtual void begin_polynomial( unsigned int variables, std::string  coeff_typename );
    virtual void end_polynomial();
    
    virtual void begin_monom( std::string coefficient );
    virtual void end_monom()
    virtual void accept_integer( std::string s );
    
    iterator objects_begin() { return object_container.begin(); }
    iterator objects_end()   { return object_container.end();   }
};

CGAL_BENCHMARK_END_NAMESPACE
