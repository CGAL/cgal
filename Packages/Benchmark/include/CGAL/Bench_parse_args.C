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
