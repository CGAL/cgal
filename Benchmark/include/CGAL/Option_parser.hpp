// Copyright (c) 1997  Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// 
// Author(s)     : Efi Fogel <efif@post.tau.ac.il>

#ifndef CGAL_BENCH_OPTION_PARSER_HPP
#define CGAL_BENCH_OPTION_PARSER_HPP

#include <boost/program_options.hpp>

#include "CGAL/benchmark_basic.hpp"

CGAL_BENCHMARK_BEGIN_NAMESPACE

namespace po = boost::program_options;

class Bench_option_parser {
public:
  /*! Constructor */
  Bench_option_parser();

  /*! Destructor */
  virtual ~Bench_option_parser();

  /*! Parse the options */
  void operator()(po::variables_map & variable_map);
  
  /*! Obtain the bench-option description */
  const po::options_description & get_opts() const { return m_bench_opts; }

  /*! Display the log-table header? */
  bool is_print_header() const { return m_print_header; }

  /*! Obtain the numbr of charactes in the name field */
  unsigned int get_name_length() const { return m_name_length; }

  /*! Obtain the number of seconds allotted to perform the bench */
  unsigned int get_seconds()  const { return m_seconds; }

  /*! Obtain the number of samples to perform the bench */
  unsigned int get_samples()  const { return m_samples; }

  /*! Obtain the number of iterations */
  unsigned int get_iterations()  const { return m_iterations; }
  
protected:
  /*! The bench options */
  po::options_description m_bench_opts;

private:
  /*! Indicates whether to display the log-table header */
  bool m_print_header;

  /*! The numbr of charactes in the name field */
  unsigned int m_name_length;

  /*! Number of seconds allotted to perform the bench */
  unsigned int m_seconds;

  /*! Number of samples to perform the bench */
  unsigned int m_samples;

  /*! Number of iterations */
  unsigned int m_iterations;  
};

CGAL_BENCHMARK_END_NAMESPACE

#endif
