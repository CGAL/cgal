#ifndef NUMBER_TYPE_H
#define NUMBER_TYPE_H

#include <CGAL/config.h>
#include "short_names.h"
#include <CGAL/basic.h>
#include "bench_config.h"

#if BENCH_NT == LEDA_REAL_NT
#include <CGAL/leda_real.h>

#elif BENCH_NT == LEDA_RAT_NT
#include <CGAL/leda_rational.h>

#elif BENCH_NT == LAZY_RATIONAL_NT
#include <CGAL/leda_rational.h>
#include <CGAL/Lazy_exact_nt.h>

#elif BENCH_NT == GMPQ_NT
#include <CGAL/Gmpq.h>

#elif BENCH_NT == LAZY_GMPQ_NT
#include <CGAL/Gmpq.h>
#include <CGAL/Lazy_exact_nt.h>

#elif BENCH_NT == QUOTIENT_MP_FLOAT_NT
#include <CGAL/MP_Float.h>
#include <CGAL/Quotient.h>

#elif BENCH_NT == QUOTIENT_GMPZ_NT
#include <CGAL/Gmpz.h>
#include <CGAL/Quotient.h>

#elif BENCH_NT == LAZY_QUOTIENT_MP_FLOAT_NT
#include <CGAL/MP_Float.h>
#include <CGAL/Quotient.h>
#include <CGAL/Lazy_exact_nt.h>

#elif BENCH_NT == DOUBLE_NT
#include "Double.h"

#else
#error No Number Type (NT) specified! 
#endif

#if BENCH_NT == LEDA_REAL_NT
typedef leda_real                                       NT;
typedef NT                                              WNT;
#define NUMBER_TYPE "Leda Real"

#elif BENCH_NT == LEDA_RAT_NT
typedef leda_rational                                   NT;
typedef NT                                              WNT;
#define NUMBER_TYPE "Leda Rat"

#elif BENCH_NT == LAZY_LEDA_RAT_NT
typedef leda_rational                                   NT;
typedef CGAL::Lazy_exact_nt<NT>                         WNT;
#define NUMBER_TYPE "Lazy Leda Rat"

#elif BENCH_NT == GMPQ_NT
typedef CGAL::Gmpq                                      NT;
typedef NT                                              WNT;
#define NUMBER_TYPE "Gmpq"

#elif BENCH_NT == LAZY_GMPQ_NT
typedef CGAL::Gmpq                                      NT;
typedef CGAL::Lazy_exact_nt<NT>                         WNT;
#define NUMBER_TYPE "Lazy Gmpq"

#elif BENCH_NT == QUOTIENT_MP_FLOAT_NT
typedef CGAL::Quotient<CGAL::MP_Float>                  NT;
typedef NT                                              WNT;
#define NUMBER_TYPE "Quotient MP Float"

#elif BENCH_NT == QUOTIENT_GMPZ_NT
typedef CGAL::Quotient<CGAL::Gmpz>                      NT;
typedef NT                                              WNT;
#define NUMBER_TYPE "Quotient Gmpz"

#elif BENCH_NT == LAZY_QUOTIENT_MP_FLOAT_NT
typedef CGAL::Quotient<CGAL::MP_Float>                  NT;
typedef CGAL::Lazy_exact_nt<NT>                         WNT;
#define NUMBER_TYPE "Lazy Quotient MP Float"

#elif BENCH_NT == DOUBLE_NT
typedef CGAL::Double                                    NT;
typedef NT                                              WNT;
#define NUMBER_TYPE "Double"

#else
#error No Number Type (NT) Specified
#endif

#endif
