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
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Efi Fogel <efif@post.tau.ac.il>
#include "Bench_parse_args.h"

CGAL_BEGIN_NAMESPACE

char Bench_parse_args::s_option_str[] = "b:d:hi:I:ps:t:v";

char * Bench_parse_args::s_io_opts[] = {"format", "f", NULL};

char * Bench_parse_args::s_format_opts[] =
  {"rat", "int", "float", "r", "i", "f", NULL};

char * Bench_parse_args::s_bench_opts[] =
  {"type_name", "type_mask", "strategy_name", "strategy_mask",
   "header", "name_length",
   "tn", "tm", "sn", "sm", "h", "nl", NULL};

char * Bench_parse_args::s_type_opts[] =
  {"increment", "aggregate", "display", "subcurves", "points",
   "i", "a", "d", "c", "p", NULL};

char * Bench_parse_args::s_strategy_opts[] =
  {"trapezoidal", "naive", "walk", "dummy", "t", "n", "w", "d"};

char * Bench_parse_args::s_bool_opts[] =
  {"true", "false", "t", "f", NULL};

CGAL_END_NAMESPACE
