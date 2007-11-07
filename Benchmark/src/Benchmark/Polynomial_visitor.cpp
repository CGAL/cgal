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
#include <algorithm>

#include <CGAL/basic.h>
#include <CGAL/assertions.h>
#include <CGAL/Polynomial.h>
#include <CGAL/Polynomial_traits_d.h>
#include <CGAL/Sqrt_extension.h>
#include <CGAL/Benchmark/Polynomial_visitor.hpp>


// work around c++ preprocessor. rationale: 
//   (1) giving a<b,c> as macro parameter directly does not work
//   (2) surrounding a<b,c> with () also does not work, because then,
//       c++ does not recognise (a<b,c>) as a type
typedef CGAL::Sqrt_extension<int,int> CGAL_Sqrt_extension;

CGAL_BENCHMARK_BEGIN_NAMESPACE


Polynomial_visitor::Polynomial_visitor() {}

void
Polynomial_visitor::begin_polynomial( unsigned int variables,
                                      std::string  coeff_typename )
{
  //std::cout << "begin polynomial with coeff type " << coeff_typename << std::endl;
  number_of_variables = variables;
  if( coeff_typename == "Integer" ) {
    coefficient_numbertype = Integer;
  } else if( coeff_typename == "Sqrt_extension<Integer,Integer>" ) {
    coefficient_numbertype = Sqrt_ext_int_int;
  } else {
    coefficient_numbertype = Not_supported;
  }
}

void Polynomial_visitor::end_polynomial() {
#define CGAL_BENCHMARK_POLYNOMIAL_VISITOR_CHECK_TYPE(NT,NAME)               \
    case NAME: {                                                               \
      CGAL::Polynomial_traits_d<CGAL::Polynomial<NT> >::Construct_polynomial construct_polynomial; \
      typedef std::pair< std::vector<int>, NT > Monom_rep;                \
      Monom_rep monom_rep;                                                     \
      std::vector< Monom_rep > polygon_rep;                                    \
      for( std::vector< Monom >::iterator it = polynomial.begin(); it != polynomial.end(); ++it ) { \
        NT coeff;                                                                \
        /*std::cout << "expecting coeff object with type " << #NAME  << std::endl;*/ \
        if( assign( coeff, it->coefficient ) )                                   \
          polygon_rep.push_back( std::make_pair( it->exponent_vector, coeff ) ); \
        else                                                                     \
          CGAL_error_msg("cannot happen");                                           \
      }                                                                          \
      object_container.push_back(                                                \
        make_object(                                                             \
          construct_polynomial( polygon_rep.begin(), polygon_rep.end() ) ) );    \
      break; }

  //std::cout << "end polynomial" << std::endl;
  switch( coefficient_numbertype ) {
    CGAL_BENCHMARK_POLYNOMIAL_VISITOR_CHECK_TYPE( int, Integer )
    CGAL_BENCHMARK_POLYNOMIAL_VISITOR_CHECK_TYPE( CGAL_Sqrt_extension, Sqrt_ext_int_int )
    default:
      CGAL_error_msg("if you can read this, something went wrong");
  }
#undef CGAL_BENCHMARK_POLYNOMIAL_VISITOR_CHECK_TYPE
  coefficient_numbertype = Not_supported;
  number_of_variables = 0;
  polynomial.clear();
}

void
Polynomial_visitor::begin_monom( std::string coefficient )
{
#define CGAL_BENCHMARK_POLYNOMIAL_VISITOR_CHECK_TYPE( NT, NAME )     \
    case NAME : {                                                    \
      NT coeff;                                                      \
      input >> coeff;                                                \
      polynomial.back().coefficient = CGAL::make_object( coeff );    \
      /*std::cout << "created coeff object with type " << #NAME  << std::endl;*/ \
      break; }
  //std::cout << "begin monom" << std::endl;
  std::stringstream input( coefficient );              
  polynomial.push_back( Polynomial_visitor::Monom() );

  switch( coefficient_numbertype ) {
    CGAL_BENCHMARK_POLYNOMIAL_VISITOR_CHECK_TYPE( int, Integer )
    CGAL_BENCHMARK_POLYNOMIAL_VISITOR_CHECK_TYPE( CGAL_Sqrt_extension, Sqrt_ext_int_int )
    default:
      CGAL_error_msg("if you can read this, something went wrong");
  }

  inside_monom = true;
#undef CGAL_BENCHMARK_POLYNOMIAL_VISITOR_CHECK_TYPE
}

void
Polynomial_visitor::end_monom()
{
  //std::cout << "end monom" << std::endl;
  std::reverse( polynomial.back().exponent_vector.begin(), polynomial.back().exponent_vector.end() );
  inside_monom = false;
}

void
Polynomial_visitor::accept_integer( std::string s )
{
  //std::cout << "accept int" << std::endl;
  if( ! inside_monom )
    return;
  polynomial.back().exponent_vector.push_back( atoi( s.c_str() ) );
}

CGAL_BENCHMARK_END_NAMESPACE
