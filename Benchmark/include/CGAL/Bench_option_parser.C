// Copyright (c) 1997  Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// 
//
// Author(s)     : Efi Fogel <efif@post.tau.ac.il>

#include <CGAL/basic.h>
#include <CGAL/Bench_option_parser.h>

CGAL_BEGIN_NAMESPACE

/*! Constructor */
Bench_option_parser::Bench_option_parser() :
  m_bench_opts("CGAL bench options"),
  m_print_header(true),
  m_name_length(32),
  m_seconds(1),
  m_samples(10),  
  m_iterations(1)
{
  // Generic options:
  m_bench_opts.add_options()
    ("header,p", po::value<bool>(&m_print_header)->default_value(true),
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
Bench_option_parser::~Bench_option_parser() {}

/*! Parse the options */
void Bench_option_parser::operator()(po::variables_map & variable_map)
{
}

CGAL_END_NAMESPACE
