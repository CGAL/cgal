#include "Bench_parse_args.h"

CGAL_BEGIN_NAMESPACE

char Bench_parse_args::s_optionStr[] = "b:d:hi:I:s:t:v";

OPTS_CONST char * Bench_parse_args::s_IOOpts[] = {"format", "f", NULL};

OPTS_CONST char * Bench_parse_args::s_formatOpts[] =
  {"rat", "int", "float", "r", "i", "f", NULL};

OPTS_CONST char * Bench_parse_args::s_benchOpts[] =
  {"type_name", "type_mask", "strategy_name", "strategy_mask",
   "header", "name_length",
   "tn", "tm", "sn", "sm", "h", "nl", NULL};

OPTS_CONST char * Bench_parse_args::s_typeOpts[] =
  {"increment", "aggregate", "display", "i", "a", "d", NULL};

OPTS_CONST char * Bench_parse_args::s_strategyOpts[] =
  {"trapezoidal", "naive", "walk", "t", "n", "w"};

OPTS_CONST char * Bench_parse_args::s_boolOpts[] =
  {"true", "false", "t", "f", NULL};

CGAL_END_NAMESPACE
