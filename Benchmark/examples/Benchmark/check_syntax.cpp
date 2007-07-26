//! \file examples/Benchmark/check_syntax.cpp
// An example that checks the syntax of an input data file

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
// Copyright (c) 2004  Max-Planck-Institut Saarbruecken (Germany)
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
// $Source: /KM/projects/ecg/CVS/BMTools/parser/src/check_syntax.C,v $
// $Revision$ $Date$
// $Name:  $
//
// Author(s) : Lutz Kettner
**************************************************************************/

#include <iostream>

#include "CGAL/Benchmark/benchmark_format.hpp"

namespace cb = CGAL::benchmark;

struct Checker : public cb::Benchmark_visitor {
    Checker() {}
    virtual void token_not_handled( std::string s) {}
    virtual void accept_benchmark_name( std::string s) {
        Benchmark_visitor::accept_benchmark_name(s);
        std::cerr << "name '" << s << "', ";
    }

    virtual void accept_classification( std::string problem,
                                        std::string geom,
                                        std::string clas,
                                        std::string family,
                                        std::string instance,
                                        std::string release) {

        if ((problem != "Arrangement") && (problem != "CSG")
            && (problem != " "))
            error_handler( "classification error");

        if ((geom != "Lines") && (geom != "Circles") && (geom != "Conics")
            && (geom != "Cubics") && (geom != "Quartics")
            && (geom != "ArbitraryDegreeCurves") && ( geom != "Quadrics")
            && (geom != "Tori") && (geom != "Planes") && (geom != " "))
            error_handler( "classification error" );

        if ((clas != "FullCurves") && (clas != "Ellipses")
            && (clas != "BoundedArcs") && (clas != "UnboundedArcs")
            && (clas != "FullSurfaces") && (clas != "BoundedPatches")
            && (clas != "UnboundedPatches") && (clas != " "))
            error_handler( "classification error" );
    }
};



int main( int argc, char* argv[] ) {
    int exit_status = 0;
    if ( argc < 2) {
        Checker check;
        std::cerr << "File read from <cin>, ";
        if ( cb::benchmark_parse_stream( std::cin, "<cin>", & check)) {
            std::cerr << "is o.k." << std::endl;
        } else {
            std::cerr << "is malformed." << std::endl;
            exit_status = 1;
        }
    } else {
        for ( int i = 1; i < argc; ++i) {
            Checker check;
            std::cerr << "File '" << argv[i] << "', ";
            if ( cb::benchmark_parse_file( argv[i], & check)) {
                std::cerr << "is o.k." << std::endl;
            } else {
                std::cerr << "is malformed." << std::endl;
                exit_status = 1;
            }
        }
    }
    return exit_status;
}

#endif
