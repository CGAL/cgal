#ifndef CGAL_TEST_NUMBER_TYPE_H
#define CGAL_TEST_NUMBER_TYPE_H

// ============================================================================
// Number type includes:
// ============================================================================
#if TEST_NT == DOUBLE_NT
#include "Double.h"

#elif TEST_NT == MP_FLOAT_NT
#include <CGAL/MP_Float.h>

#elif TEST_NT == GMPZ_NT
#include <CGAL/gmpxx.h>

#elif TEST_NT == LEDA_RAT_NT
#include <CGAL/leda_rational.h>

#elif TEST_NT == LAZY_LEDA_RAT_NT
#include <CGAL/leda_rational.h>
#include <CGAL/Lazy_exact_nt.h>

#elif TEST_NT == GMPQ_NT
#include <gmpxx.h>

#elif TEST_NT == CGAL_GMPQ_NT
#include <CGAL/Gmpq.h>
#include <CGAL/Gmpz.h>

#elif TEST_NT == LAZY_CGAL_GMPQ_NT
#include <CGAL/Gmpq.h>
#include <CGAL/Gmpz.h>
#include <CGAL/Lazy_exact_nt.h>

#elif TEST_NT == QUOTIENT_MP_FLOAT_NT
#include <CGAL/MP_Float.h>
#include <CGAL/Quotient.h>

#elif TEST_NT == QUOTIENT_CGAL_GMPZ_NT
#include <CGAL/Gmpz.h>
#include <CGAL/Quotient.h>

#elif TEST_NT == LAZY_QUOTIENT_MP_FLOAT_NT
#include <CGAL/MP_Float.h>
#include <CGAL/Quotient.h>
#include <CGAL/Lazy_exact_nt.h>

#elif TEST_NT == LEDA_REAL_NT
#include <CGAL/leda_real.h>

#elif TEST_NT == CORE_EXPR_NT
#include <CGAL/CORE_Expr.h>

#elif TEST_NT == LAZY_GMPZ_NT
#include <CGAL/Lazy_exact_nt.h>
#include <CGAL/gmpxx.h>

#elif TEST_NT == LEDA_INT_NT
#include <CGAL/leda_integer.h>

#elif TEST_NT == CGAL_GMPZ_NT
#include <CGAL/Gmpz.h>

#elif TEST_NT == CORE_INT_NT
#include <CGAL/CORE_BigInt.h>

#elif TEST_NT == CORE_RAT_NT          //OS - new
#include <CGAL/Arithmetic_kernel.h>

#else
#error No Number Type (NT) specified!
#endif

// ============================================================================
// Number type typedef:
// ============================================================================
#if TEST_NT == DOUBLE_NT
typedef CGAL::Double                                    Basic_number_type;
typedef Basic_number_type                               Number_type;
typedef Basic_number_type                               Ring_type;
#define NUMBER_TYPE "Double"

#elif TEST_NT == MP_FLOAT_NT
typedef CGAL::MP_Float                                  Basic_number_type;
typedef Basic_number_type                               Number_type;
typedef Basic_number_type                               Ring_type;
#define NUMBER_TYPE "MP Float"

#elif TEST_NT == GMPZ_NT
typedef ::mpz_class                                     Basic_number_type;
typedef Basic_number_type                               Number_type;
typedef Basic_number_type                               Ring_type;
#define NUMBER_TYPE "Gmpz"

#elif TEST_NT == LEDA_RAT_NT
typedef leda_rational                                   Basic_number_type;
typedef Basic_number_type                               Number_type;
typedef leda_integer                                    Ring_type;
#define NUMBER_TYPE "Leda Rat"

#elif TEST_NT == LAZY_LEDA_RAT_NT
typedef leda_rational                                   Basic_number_type;
typedef CGAL::Lazy_exact_nt<NT>                         Number_type;
typedef leda_integer                                    Ring_type;
#define NUMBER_TYPE "Lazy Leda Rat"

#elif TEST_NT == GMPQ_NT
typedef ::mpq_class                                     Basic_number_type;
typedef Basic_number_type                               Number_type;
typedef ::mpz_class                                     Ring_type;
typedef ::mpq_class                                     FT;
#define NUMBER_TYPE "Gmpq"

#elif TEST_NT == CGAL_GMPQ_NT
typedef CGAL::Gmpq                                      Basic_number_type;
typedef Basic_number_type                               Number_type;
typedef CGAL::Gmpz                                      Ring_type;
typedef CGAL::Gmpq                                      FT;
#define NUMBER_TYPE "Cgal Gmpq"

#elif TEST_NT == LAZY_CGAL_GMPQ_NT
typedef CGAL::Gmpq                                      Basic_number_type;
typedef CGAL::Lazy_exact_nt<NT>                         Number_type;
typedef CGAL::Gmpz                                      Ring_type;
typedef CGAL::Gmpq                                      FT;
#define NUMBER_TYPE "Lazy Cgal Gmpq"

#elif TEST_NT == QUOTIENT_MP_FLOAT_NT
typedef CGAL::Quotient<CGAL::MP_Float>                  Basic_number_type;
typedef Basic_number_type                               Number_type;
typedef CGAL::MP_Float                                  Ring_type;
#define NUMBER_TYPE "Quotient MP Float"

#elif TEST_NT == QUOTIENT_CGAL_GMPZ_NT
typedef CGAL::Quotient<CGAL::Gmpz>                      Basic_number_type;
typedef Basic_number_type                               Number_type;
typedef CGAL::Gmpz                                      Ring_type;
#define NUMBER_TYPE "Quotient Gmpz"

#elif TEST_NT == LAZY_QUOTIENT_MP_FLOAT_NT
typedef CGAL::Quotient<CGAL::MP_Float>                  Basic_number_type;
typedef CGAL::Lazy_exact_nt<NT>                         Number_type;
typedef CGAL::MP_Float                                  Ring_type;
#define NUMBER_TYPE "Lazy Quotient MP Float"

#elif TEST_NT == LEDA_REAL_NT
typedef leda_real                                       Basic_number_type;
typedef Basic_number_type                               Number_type;
typedef Basic_number_type                               Ring_type;
#define NUMBER_TYPE "Leda Real"

#elif TEST_NT == CORE_EXPR_NT
typedef CORE::Expr                                      Basic_number_type;
typedef Basic_number_type                               Number_type;
typedef Basic_number_type                               Ring_type;
#define NUMBER_TYPE "Core Expr"

#elif TEST_NT == LAZY_GMPZ_NT
typedef ::mpz_class                                     Basic_number_type;
typedef CGAL::Lazy_exact_nt<NT>                         Number_type;
typedef ::mpz_class                                     Ring_type;
#define NUMBER_TYPE "Lazy Gmpz"

#elif TEST_NT == LEDA_INT_NT
typedef leda_integer                                    Basic_number_type;
typedef Basic_number_type                               Number_type;
typedef Basic_number_type                               Ring_type;
#define NUMBER_TYPE "Leda Int"

#elif TEST_NT == CGAL_GMPZ_NT
typedef CGAL::Gmpz                                      Basic_number_type;
typedef Basic_number_type                               Number_type;
typedef Basic_number_type                               Ring_type;
#define NUMBER_TYPE "Cgal Gmpz"

#elif TEST_NT == CORE_INT_NT
typedef CORE::BigInt                                    Basic_number_type;
typedef Basic_number_type                               Number_type;
typedef Basic_number_type                               Ring_type;
#define NUMBER_TYPE "CORE BigInt"

#elif TEST_NT == CORE_RAT_NT
typedef CGAL::CORE_arithmetic_kernel::Rational                    Basic_number_type;
typedef Basic_number_type                               Number_type;
typedef Basic_number_type                               Ring_type;
#define NUMBER_TYPE "CORE Rational"

#else
#error No Number Type (NT) Specified
#endif

#endif
