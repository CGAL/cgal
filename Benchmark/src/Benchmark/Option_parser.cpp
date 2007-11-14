// Copyright (c) 1997  Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it and/or
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
// Author(s)     : Efi Fogel <efif@post.tau.ac.il>

#include "CGAL/Benchmark/config.hpp"
#include "CGAL/Benchmark/Option_parser.hpp"

CGAL_BENCHMARK_BEGIN_NAMESPACE

/*! Constructor */
 Option_parser:: Option_parser() :
  m_bench_opts("CGAL bench options"),
  m_print_header(true),
  m_name_length(32),
  m_seconds(1),
  m_samples(10),  
  m_iterations(1)
{
  // Generic options:
  m_bench_opts.add_options()
    ("print-header,p", po::value<bool>(&m_print_header)->default_value(true),
     "print header")
    ("name-length,n",
     po::value<unsigned int>(&m_name_length)->default_value(32),
     "name-field length")
    ("samples,s", po::value<unsigned int>(&m_samples)->default_value(10),
     "number of samples")
    ("seconds,t", po::value<unsigned int>(&m_seconds)->default_value(1),
     "number of seconds")
    ("iterations,i", po::value<unsigned int>(&m_iterations)->default_value(1),
     "number of iterations")
    ;
}

/*! Destructor */
 Option_parser::~ Option_parser() {}

/*! Parse the options */
void  Option_parser::operator()(po::variables_map & variable_map)
{
}

CGAL_BENCHMARK_END_NAMESPACE
