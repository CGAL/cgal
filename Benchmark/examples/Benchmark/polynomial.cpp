//! \file examples/Benchmark/polynomial.cpp
// Measure the performance of ???

#include <CGAL/basic.h>

#ifndef CGAL_USE_BOOST_PROGRAM_OPTIONS
#include <iostream>
int main()
{
  std::cout << "Sorry, this example needs boost program options ..."
            << std::endl;
  return 0;
}
#else

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
// $Source: $
// $Revision: 34400 $ $Date: 2006-09-20 12:52:41 +0200 (Wed, 20 Sep 2006) $
// $Name:  $
//
// Author(s) : Andreas Meyer
**************************************************************************/

#include <iostream>
#include <string>
#include <sstream>

#include <CGAL/Benchmark/benchmark_format.hpp>
#include <CGAL/Sqrt_extension.h>

namespace cb = CGAL::benchmark;



struct Polynomial_reader: public cb::Benchmark_visitor {
    unsigned char current_variable_name;
    bool inside_monom;
    bool first_monom;
    enum Coefficient_numbertype { Not_supported = 0, Integer, Sqrt_ext_int_int };
    
    Coefficient_numbertype  coefficient_numbertype;
    

    Polynomial_reader() : inside_monom( false ) {}

    virtual void
    token_not_handled( std::string s)
    {}

    // accept everything
    virtual void accept_classification( std::string problem,
                                        std::string geom,
                                        std::string clas,
                                        std::string family,
                                        std::string instance,
                                        std::string release) {}


    virtual void
    begin_polynomial( unsigned int variables,
                      std::string  coeff_typename )
    {
        std::cout << "number of variables: " << variables << std::endl;
        std::cout << "typename: " << coeff_typename << std::endl;
        if( coeff_typename == "Sqrt_extension<Integer,Integer>" ) {
          coefficient_numbertype = Sqrt_ext_int_int;
        } else if( coeff_typename == "Integer" ) {
          coefficient_numbertype = Integer;
        } else {
          std::cout << "coefficient numbertype " << coeff_typename << " not supported." << std::endl;
          coefficient_numbertype = Not_supported;
        }
        first_monom = true;
    }

    virtual void
    end_polynomial()
    {
      std::cout << std::endl;
    }

    virtual void
    begin_monom( std::string coefficient )
    {
        switch( coefficient_numbertype ) {
          case Sqrt_ext_int_int: {
            //std::cout << "parsing coefficient " << coefficient << std::endl;
            CGAL::Sqrt_extension<int,int> ex;
            std::stringstream input( coefficient );
            input >> ex;
            break;
          }
          case Integer:
            break;
          default: 
            ;
        }
        if( ! first_monom ) {
           std::cout << " + ";
        }
        first_monom = false;
        if( coefficient == "1" )
          coefficient = "";
        else if( coefficient == "-1" )
          coefficient = "-";
        std::cout << coefficient << "(";
        current_variable_name = 'a';
        inside_monom = true;
    }

    virtual void
    end_monom()
    {
        std::cout << ")";
        inside_monom = false;
    }

    virtual void
    accept_integer( std::string s )
    {
        if( ! inside_monom )
            return;
        if( s != "0" )
          std::cout << current_variable_name << "^" << s;
        ++current_variable_name;
    }
};

int main( int argc, char* argv[] ) {
    int exit_status = 0;
    if ( argc < 2) {
        Polynomial_reader reader;
        if ( !cb::benchmark_parse_stream( std::cin, "<cin>", & reader)) {
            std::cerr << "input malformed." << std::endl;
            exit_status = 1;
        }
    }
    return exit_status;
}

#endif
