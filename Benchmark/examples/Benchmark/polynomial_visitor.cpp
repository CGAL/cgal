//! \file examples/Benchmark/polynomial_visitor.cpp
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
// $Revision$ $Date$
// $Name:  $
//
// Author(s) : Andreas Meyer
**************************************************************************/

#include <iostream>
#include <string>
#include <sstream>

#include <CGAL/Benchmark/benchmark_format.hpp>
#include <CGAL/Benchmark/Polynomial_visitor.hpp>

#include <CGAL/Sqrt_extension.h>
#include <CGAL/Polynomial.h>

namespace cb = ::CGAL::benchmark;


int main( int argc, char* argv[] ) {
    int exit_status = 0;
    if ( argc < 2) {
        cb::Polynomial_visitor pv;
        if ( cb::benchmark_parse_stream( std::cin, "<cin>", & pv)) {
	    std::cout << "ok" << std::endl;
	    typedef cb::Polynomial_visitor::iterator iterator;
            CGAL::Polynomial< int                           > p_int;
            CGAL::Polynomial< CGAL::Sqrt_extension<int,int> > p_sqrt_ext;
	    for( iterator it = pv.objects_begin(); it != pv.objects_end(); ++it ) {
	      if( CGAL::assign( p_sqrt_ext, *it ) ) {
	        std::cout << "read CGAL::Polynomial< CGAL::Sqrt_extension<int,int> >: " << p_sqrt_ext << std::endl;
	      }
	      else if( CGAL::assign( p_int, *it ) ) {
	        std::cout << "read CGAL::Polynomial< int >: " << p_int << std::endl;
	      } else {
	        std::cout << "read unknown object" << std::endl;
	      }
	    }
	    
	} else {
            std::cerr << "input malformed." << std::endl;
            exit_status = 1;
        }
    }
    return exit_status;
}

#endif
