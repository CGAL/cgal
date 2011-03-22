#ifndef NUMBER_TYPE_HPP
#define NUMBER_TYPE_HPP

#include <CGAL/config.h>
#include <CGAL/basic.h>

#include "bench_config.hpp"

#if BENCH_NT == DOUBLE_NT
#include "Double.hpp"

#elif BENCH_NT == MP_FLOAT_NT
#include <CGAL/MP_Float.h>

#elif BENCH_NT == GMPZ_NT
#include <CGAL/gmpxx.h>

#elif BENCH_NT == LEDA_RAT_NT
#include <CGAL/leda_rational.h>

#elif BENCH_NT == LAZY_LEDA_RAT_NT
#include <CGAL/leda_rational.h>
#include <CGAL/Lazy_exact_nt.h>

#elif BENCH_NT == GMPQ_NT
#include <gmpxx.h>

#elif BENCH_NT == CGAL_GMPQ_NT
#include <CGAL/Gmpq.h>
#include <CGAL/Gmpz.h>

#elif BENCH_NT == LAZY_CGAL_GMPQ_NT
#include <CGAL/Gmpq.h>
#include <CGAL/Gmpz.h>
#include <CGAL/Lazy_exact_nt.h>

#elif BENCH_NT == QUOTIENT_MP_FLOAT_NT
#include <CGAL/MP_Float.h>
#include <CGAL/Quotient.h>

#elif BENCH_NT == QUOTIENT_CGAL_GMPZ_NT
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
#if ! defined(CGAL_USE_LEDA)
#error "Leda not supported"
#endif

#elif BENCH_NT == NIX_CORE_FIELD_WITH_SQRT_NT
#include <NiX/Arithmetic_traits.h>
#if !defined(CGAL_USE_CORE)
#error "Core not supported!"
#endif

#elif BENCH_NT == LAZY_GMPZ_NT
#include <CGAL/Lazy_exact_nt.h>
#include <CGAL/gmpxx.h>

#else
#error No Number Type (NT) specified! 
#endif

// Typedefs:

#if BENCH_NT == DOUBLE_NT
typedef CGAL::Double                                          Number_type;
#define NUMBER_TYPE "Double"

#elif BENCH_NT == MP_FLOAT_NT
typedef CGAL::MP_Float                                        Number_type;
#define NUMBER_TYPE "MP Float"

#elif BENCH_NT == GMPZ_NT
typedef ::mpz_class                                           Number_type;
#define NUMBER_TYPE "Gmpz"

#elif BENCH_NT == LEDA_RAT_NT
typedef leda_rational                                         Number_type;
#define NUMBER_TYPE "Leda Rat"

#elif BENCH_NT == LAZY_LEDA_RAT_NT
typedef CGAL::Lazy_exact_nt<leda_rational>                    Number_type;
#define NUMBER_TYPE "Lazy Leda Rat"

#elif BENCH_NT == GMPQ_NT
typedef ::mpq_class                                           Number_type;
#define NUMBER_TYPE "Gmpq"

#elif BENCH_NT == CGAL_GMPQ_NT
typedef CGAL::Gmpq                                            Number_type;
#define NUMBER_TYPE "Cgal Gmpq"

#elif BENCH_NT == LAZY_CGAL_GMPQ_NT
typedef CGAL::Lazy_exact_nt<CGAL::Gmpq>                       Number_type;
#define NUMBER_TYPE "Lazy Cgal Gmpq"

#elif BENCH_NT == QUOTIENT_MP_FLOAT_NT
typedef CGAL::Quotient<CGAL::MP_Float>                        Number_type;
#define NUMBER_TYPE "Quotient MP Float"

#elif BENCH_NT == QUOTIENT_CGAL_GMPZ_NT
typedef CGAL::Quotient<CGAL::Gmpz>                            Number_type;
#define NUMBER_TYPE "Quotient Gmpz"

#elif BENCH_NT == LAZY_QUOTIENT_MP_FLOAT_NT
typedef CGAL::Lazy_exact_nt<CGAL::Quotient<CGAL::MP_Float> >  Number_type;
#define NUMBER_TYPE "Lazy Quotient MP Float"

#elif BENCH_NT == LEDA_REAL_NT
typedef leda_real                                             Number_type;
#define NUMBER_TYPE "Leda Real"

#elif BENCH_NT == CORE_EXPR_NT
typedef CORE::Expr                                            Number_type;
#define NUMBER_TYPE "Core Expr"

#elif BENCH_NT == NIX_LEDA_FIELD_WITH_SQRT_NT
typedef NiX::LEDA_arithmetic_traits                           Arithmetic_traits;
typedef Arithmetic_traits::Field_with_sqrt                    Number_type;
#define NUMBER_TYPE "NiX Leda Real"

#elif BENCH_NT == NIX_CORE_FIELD_WITH_SQRT_NT
typedef NiX::CORE_arithmetic_traits                           Arithmetic_traits;
typedef Arithmetic_traits::Field_with_sqrt                    Number_type;
#define NUMBER_TYPE "NiX Core Expr"

#elif BENCH_NT == LAZY_GMPZ_NT
typedef CGAL::Lazy_exact_nt<::mpz_class>                      Number_type;
#define NUMBER_TYPE "Lazy Gmpz"

#else
#error No Number Type (NT) Specified
#endif

#endif
