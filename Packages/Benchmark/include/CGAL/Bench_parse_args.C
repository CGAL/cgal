// ======================================================================
//
// Copyright (c) 1997 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       : $$
// release_date  : $$
//
// file          : include/CGAL/Bench_parse_args.C
// package       : Planar_map (5.87)
// maintainer    : Efi Fogel <efif@post.tau.ac.il>
// source        : 
// revision      : 
// revision_date : 
// author(s)     : Efi Fogel <efif@post.tau.ac.il>
//
// coordinator   : Tel-Aviv University (Dan Halperin <halperin@math.tau.ac.il>)
//
// Chapter       : 
// ======================================================================
#include "Bench_parse_args.h"

CGAL_BEGIN_NAMESPACE

char Bench_parse_args::s_option_str[] = "b:d:hi:I:s:t:v";

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
