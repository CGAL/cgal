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

#include <iostream>

#include <CGAL/Polynomial.h>
#include <CGAL/Sqrt_extension.h>
#include <CGAL/Benchmark/Polynomial_visitor.hpp>

CGAL_BENCHMARK_BEGIN_NAMESPACE

template< typename NT, int vars >
struct Polynomial_type_creator {
  typedef typename Polynomial_type_creator<NT,vars-1>::result result;
};

template< typename NT >
struct Polynomial_type_creator<NT,1> {
  typedef CGAL::Polynomial< NT > result;
};

template< typename NT >
struct Polynomial_type_creator<NT,0> {
  typedef void result;
};

Polynomial_visitor::Polynomial_visitor() {

}

void
Polynomial_visitor::begin_polynomial( unsigned int variables,
                                      std::string  coeff_typename )
{
  number_of_variables = variables;
  if( coeff_typename == "Integer" ) {
    coefficient_numbertype = Integer;
    CGAL::Polynomial< int > p;
    object_container.push_back(  make_object( p ) );
  } else if( coeff_typename == "Sqrt_extension<Integer,Integer>" ) {
    coefficient_numbertype = Sqrt_ext_int_int;
    CGAL::Polynomial< CGAL::Sqrt_extension<int,int> > p;
    object_container.push_back(  make_object( p ) );
  } else {
    coefficient_numbertype = Not_supported;
  }
}

void Polynomial_visitor::end_polynomial() {}

void
Polynomial_visitor::begin_monom( std::string coefficient )
{

#define CGAL_BENCHMARK_POLYNOMIAL_VISITOR_CHECK_TYPE( NT )             \
    if( CGAL::Polynomial< NT > *p =                                    \
        object_cast< CGAL::Polynomial< NT > >( object_container.back() ) ) { \
      NT c;                                                            \
      std::stringstream input( coefficient );                          \
      input >> c;                                                      \
      coefficient_value = CGAL::make_object( ex );                     \
    }

    exponent_vector.clear();
    
    CGAL_BENCHMARK_POLYNOMIAL_VISITOR_CHECK_TYPE( int )
    CGAL_BENCHMARK_POLYNOMIAL_VISITOR_CHECK_TYPE( CGAL::Sqrt_extension<int,int> )
    // ...
    
    inside_monom = true;
#undef CGAL_BENCHMARK_POLYNOMIAL_VISITOR_CHECK_TYPE
}

void
Polynomial_visitor::end_monom()
{
    exponent_vector.reverse();
    inside_monom = false;
}

void
Polynomial_visitor::accept_integer( std::string s )
{
    if( ! inside_monom )
        return;
    exponent_vector.push_back( atoi( s.c_str() ) );
}
