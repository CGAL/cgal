// Copyright (c) 2004  Max-Planck-Institut Saarbruecken (Germany)
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
//
// Author(s) : Lutz Kettner

#ifndef CGAL_BENCHMARK_BENCHMARK_FORMAT_HPP
#define CGAL_BENCHMARK_BENCHMARK_FORMAT_HPP

#include <string>

#include "CGAL/Benchmark/config.hpp"
#include "CGAL/Benchmark/Benchmark_visitor.hpp"

CGAL_BENCHMARK_BEGIN_NAMESPACE

// Magic name and file version number
// ----------------------------------
// The name and the version number <major.minor> of the benchmark file format.
// If the file has a different format name, it is rejected. If the
// Major number in the file is higher than the major number listed
// here then the file is also rejected. If the <major.minor> number
// in the file is lower than the the number listed here, the parser
// assumes backwards compatibility and continues parsing.
const std::string BENCHMARK_FORMAT_NAME( "AcsBenchmark");
const int         BENCHMARK_FORMAT_MAJOR = 0;
const int         BENCHMARK_FORMAT_MINOR = 1;


// Public interface of the parser component
// ----------------------------------------

// Opens file 'name' and parses it. Uses visitor 'v' while parsing.
// Returns false if something went wrong. See the visitor for details
// of the error reporting.
bool benchmark_parse_file( std::string name, Benchmark_visitor* v);

// Starts parsing from stream 'in' with the associated filename 'name'
// (or analogous meaning for different streams) counting linenumbers
// starting from 'n'. Uses visitor 'v' while parsing. Returns false if
// something went wrong. See the visitor for details of the error reporting.
bool benchmark_parse_stream( std::istream&  in, std::string name, 
                             Benchmark_visitor* v, int n = 1);


// Public interface of the lexer component
// ---------------------------------------
// Current input stream and associated values. READ ONLY! Do not modify them!
// The lexer works with internal buffers and would not sync properly!
extern int           benchmark_linenumber;
extern std::string   benchmark_filename;

// Initialize lexer to scan from input stream with name and linenumber.
// The caller is responsible for the lifetime of 'in' that must live
// during the scan. Used by benchmark_parse_file and does not have to 
// be called separately.
void benchmark_init_lexer( 
    std::istream& in, std::string name, int linenumber = 1);

// Writes a trace of the current include file nesting to the 'out' stream.
// Appends the 'fill' string after each file listed.
void benchmark_include_file_trace( std::ostream& out, std::string fill);

CGAL_BENCHMARK_END_NAMESPACE

#endif // BENCHMARK_PARSER_H //
