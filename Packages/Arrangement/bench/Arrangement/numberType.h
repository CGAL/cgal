#ifndef NUMBER_TYPE_H
#define NUMBER_TYPE_H

#include <CGAL/config.h>
#include "short_names.h"
#include <CGAL/basic.h>
#include "bench_config.h"

#if BENCH_NT == DOUBLE_NT
#include "Double.h"

#elif BENCH_NT == GMPZ_NT
#include <CGAL/gmpxx.h>

#elif BENCH_NT == LEDA_RAT_NT
#include <CGAL/leda_rational.h>

#elif BENCH_NT == LAZY_LEDA_RAT_NT
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

#elif BENCH_NT == LEDA_REAL_NT
#include <CGAL/leda_real.h>

#elif BENCH_NT == CORE_EXPR_NT
#include <CGAL/CORE_Expr.h>

// #elif BENCH_NT == NIX_ARITHMETIC_TRAITS_NT
// #include <LiS/file_io.h>
// #include <NiX/Arithmetic_traits.h>

#elif BENCH_NT == NIX_LEDA_FIELD_WITH_SQRT_NT
#include <NiX/Arithmetic_traits.h>
#if !LiS_HAVE_LEDA
#error "Leda not supported"
#endif

#elif BENCH_NT == NIX_CORE_FIELD_WITH_SQRT_NT
#include <NiX/Arithmetic_traits.h>
#if !LiS_HAVE_CORE
#error "Core not supported!"
#endif

#else
#error No Number Type (NT) specified! 
#endif

// Typedefs:

#if BENCH_NT == DOUBLE_NT
typedef CGAL::Double                                    NT;
typedef NT                                              WNT;
#define NUMBER_TYPE "Double"

#elif BENCH_NT == GMPZ_NT
typedef ::mpz_class                                     NT;
typedef NT                                              WNT;
typedef ::mpz_class                                     RT;
typedef ::mpq_class                                     FT;
#define NUMBER_TYPE "Gmpz"

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

#elif BENCH_NT == LEDA_REAL_NT
typedef leda_real                                       NT;
typedef NT                                              WNT;
#define NUMBER_TYPE "Leda Real"

#elif BENCH_NT == CORE_EXPR_NT
typedef CORE::Expr                                      NT;
typedef NT                                              WNT;
#define NUMBER_TYPE "Core Expr"

// #elif BENCH_NT == NIX_ARITHMETIC_TRAITS_NT
// typedef NiX::Arithmetic_traits                       Arithmetic_traits;
// #if LiS_HAVE_LEDA
// typedef NiX::LEDA_arithmetic_traits Arithmetic_traits;
// #define NUMBER_TYPE "LEDA_arithmetic_traits"
// #else
// #if LiS_HAVE_CORE
// typedef NiX::CORE_arithmetic_traits Arithmetic_traits;
// #define NUMBER_TYPE "CORE_arithmetic_traits"
// #endif
// typedef Arithmetic_traits::Field_with_sqrt              NT;
// typedef NT                                              WNT;

#elif BENCH_NT == NIX_LEDA_FIELD_WITH_SQRT_NT
typedef NiX::LEDA_arithmetic_traits                     Arithmetic_traits;
typedef Arithmetic_traits::Field_with_sqrt              NT;
typedef NT                                              WNT;
#define NUMBER_TYPE "NiX Leda Real"

#elif BENCH_NT == NIX_CORE_FIELD_WITH_SQRT_NT
typedef NiX::CORE_arithmetic_traits                     Arithmetic_traits;
typedef Arithmetic_traits::Field_with_sqrt              NT;
typedef NT                                              WNT;
#define NUMBER_TYPE "NiX Core Expr"

#else
#error No Number Type (NT) Specified
#endif

#endif
