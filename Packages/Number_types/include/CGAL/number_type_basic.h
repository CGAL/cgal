// Copyright (c) 1999  Utrecht University (The Netherlands),
// ETH Zurich (Switzerland), Freie Universitaet Berlin (Germany),
// INRIA Sophia-Antipolis (France), Martin-Luther-University Halle-Wittenberg
// (Germany), Max-Planck-Institute Saarbruecken (Germany), RISC Linz (Austria),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
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
// Author(s)     : Stefan Schirra
 

#ifndef CGAL_NUMBER_TYPE_BASIC_H
#define CGAL_NUMBER_TYPE_BASIC_H

#define CGAL_PI 3.14159265358979323846

#define CGAL_NTS CGAL::
// #define CGAL_NTS CGAL::NTS::

#if ((__GNUC__ == 2) && (__GNUC_MINOR__ == 95))
#include <cmath>
#endif

// CGAL uses std::min and std::max

#include <algorithm>

CGAL_BEGIN_NAMESPACE

using std::min;
using std::max;

CGAL_END_NAMESPACE

#include <CGAL/Number_type_traits.h>
#include <CGAL/number_utils.h>
#include <CGAL/double.h>
#include <CGAL/float.h>
#include <CGAL/int.h>

// Including all number type files is necessary for compilers implementing
// two-stage name lookup (like g++ >= 3.4).
// A nicer solution needs more thought.

#ifdef CGAL_CFG_HAS_TWO_STAGE_NAME_LOOKUP

CGAL_BEGIN_NAMESPACE

// MP_Float

class MP_Float;

Comparison_result compare(const MP_Float&, const MP_Float&);
Sign sign(const MP_Float &);
MP_Float square(const MP_Float&);
MP_Float sqrt(const MP_Float&);
double to_double(const MP_Float&);
std::pair<double,double> to_interval(const MP_Float &);
bool is_finite(const MP_Float &);
bool is_valid(const MP_Float &);

// long

double to_double(long);
std::pair<double,double> to_interval(const long &);
bool is_finite(long);
bool is_valid(long);

// long long

#ifdef CGAL_USE_LONG_LONG
double to_double(long long);
std::pair<double,double> to_interval(const long long &);
bool is_finite(long long);
bool is_valid(long long);
long long int div(long long int, long long int);
#endif // CGAL_USE_LONG_LONG

// Fixed_precision_nt

class Fixed_precision_nt;

double to_double(Fixed_precision_nt);
bool is_finite(Fixed_precision_nt);
bool is_valid(Fixed_precision_nt);
std::pair<double,double> to_interval(Fixed_precision_nt);

// Quotient

template <typename> class Quotient;

template <class NT>
Quotient<NT> sqrt(const Quotient<NT> &);

template <class NT>
Comparison_result compare(const Quotient<NT>&, const Quotient<NT>&);

template <class NT>
double to_double(const Quotient<NT>&);

template <class NT>
std::pair<double,double> to_interval (const Quotient<NT>&);

template <class NT>
bool is_valid(const Quotient<NT>&);

template <class NT>
bool is_finite(const Quotient<NT>&);

// Lazy_exact_nt

template <typename> class Lazy_exact_nt;

template <typename ET>
double to_double(const Lazy_exact_nt<ET> &);

template <typename ET>
std::pair<double,double> to_interval(const Lazy_exact_nt<ET> &);

template <typename ET>
Sign sign(const Lazy_exact_nt<ET> &);

template <typename ET>
Comparison_result
compare(const Lazy_exact_nt<ET> &, const Lazy_exact_nt<ET> &);

template <typename ET>
Lazy_exact_nt<ET> abs(const Lazy_exact_nt<ET> &);

template <typename ET>
Lazy_exact_nt<ET> square(const Lazy_exact_nt<ET> &);

template <typename ET>
Lazy_exact_nt<ET> sqrt(const Lazy_exact_nt<ET> &);

template <typename ET>
Lazy_exact_nt<ET> min(const Lazy_exact_nt<ET> &, const Lazy_exact_nt<ET> &);

template <typename ET>
Lazy_exact_nt<ET> max(const Lazy_exact_nt<ET> &, const Lazy_exact_nt<ET> &);

template <typename ET>
bool is_finite(const Lazy_exact_nt<ET> &);

template <typename ET>
bool is_valid(const Lazy_exact_nt<ET> &);

// Interval_nt

template <bool> class Interval_nt;

template <bool Protected>
double to_double (const Interval_nt<Protected> &);

template <bool Protected>
std::pair<double, double> to_interval (const Interval_nt<Protected> &);

template <bool Protected>
bool is_valid (const Interval_nt<Protected> &);

template <bool Protected>
bool is_finite (const Interval_nt<Protected> &);

template <bool Protected>
Interval_nt<Protected> sqrt (const Interval_nt<Protected> &);

template <bool Protected>
Interval_nt<Protected>
min (const Interval_nt<Protected> &, const Interval_nt<Protected> &);

template <bool Protected>
Interval_nt<Protected>
max (const Interval_nt<Protected> &, const Interval_nt<Protected> &);

template <bool Protected>
Interval_nt<Protected> square (const Interval_nt<Protected> &);

template <bool Protected>
Interval_nt<Protected> abs (const Interval_nt<Protected> &);

template <bool Protected>
Sign sign (const Interval_nt<Protected> &);

template <bool Protected>
Comparison_result
compare (const Interval_nt<Protected> &, const Interval_nt<Protected> &);

// Filtered_exact

template < class, class, class, bool, class > class Filtered_exact;
struct Dynamic;

#ifndef CGAL_DENY_INEXACT_OPERATIONS_ON_FILTER
template < class CT, class ET, bool Protected, class Cache >
Filtered_exact<CT,ET,Dynamic,Protected,Cache>
div (const Filtered_exact<CT,ET,Dynamic, Protected,Cache>&,
     const Filtered_exact<CT,ET, Dynamic, Protected,Cache>&);

template < class CT, class ET, class Type, bool Protected, class Cache >
Filtered_exact<CT, ET, Type, Protected, Cache>
sqrt (const Filtered_exact<CT, ET, Type, Protected, Cache>&);

template < class CT, class ET, bool Protected, class Cache >
Filtered_exact<CT,ET,Dynamic,Protected,Cache>
gcd (const Filtered_exact<CT,ET,Dynamic, Protected,Cache>&,
     const Filtered_exact<CT,ET,Dynamic, Protected,Cache>&);

template < class CT, class ET, bool Protected, class Cache >
Filtered_exact<CT,ET,Dynamic,Protected,Cache>
square (const Filtered_exact<CT, ET, Dynamic, Protected, Cache>&);

#endif // CGAL_DENY_INEXACT_OPERATIONS_ON_FILTER

template < class CT, class ET, class Type, bool Protected, class Cache >
bool is_valid (const Filtered_exact<CT, ET, Type, Protected, Cache>&);

template < class CT, class ET, class Type, bool Protected, class Cache >
bool is_finite (const Filtered_exact<CT, ET, Type, Protected, Cache>&);

template < class CT, class ET, class Type, bool Protected, class Cache >
double to_double (const Filtered_exact<CT, ET, Type, Protected, Cache>&);

template < class CT, class ET, class Type, bool Protected, class Cache >
std::pair<double, double>
to_interval (const Filtered_exact<CT, ET, Type, Protected, Cache>&);

template < class CT, class ET, bool Protected, class Cache >
Sign sign (const Filtered_exact<CT, ET, Dynamic, Protected, Cache>&);

template < class CT, class ET, bool Protected, class Cache >
Comparison_result
compare (const Filtered_exact<CT, ET, Dynamic, Protected, Cache>&,
         const Filtered_exact<CT, ET, Dynamic, Protected, Cache>&);

template < class CT, class ET, bool Protected, class Cache >
Filtered_exact<CT,ET,Dynamic,Protected,Cache>
abs (const Filtered_exact<CT, ET, Dynamic, Protected, Cache>&);

template < class CT, class ET, bool Protected, class Cache >
Filtered_exact<CT,ET,Dynamic,Protected,Cache>
min (const Filtered_exact<CT, ET, Dynamic, Protected, Cache>&,
     const Filtered_exact<CT, ET, Dynamic, Protected, Cache>&);

template < class CT, class ET, bool Protected, class Cache >
Filtered_exact<CT,ET,Dynamic,Protected,Cache>
max (const Filtered_exact<CT, ET, Dynamic, Protected, Cache>&,
     const Filtered_exact<CT, ET, Dynamic, Protected, Cache>&);

#ifdef CGAL_USE_GMP
// Gmpz

class Gmpz;

double to_double(const Gmpz&);
Sign sign(const Gmpz &);
bool is_valid(const Gmpz &);
bool is_finite(const Gmpz &);
Gmpz sqrt(const Gmpz &);
Gmpz div(const Gmpz &, const Gmpz &);
Gmpz gcd(const Gmpz &, const Gmpz &);
Gmpz gcd(const Gmpz &, int);
std::pair<double, double> to_interval (const Gmpz &);

// Gmpq

class Gmpq;

double to_double(const Gmpq&);
Sign sign(const Gmpq &);
bool is_valid(const Gmpq &);
bool is_finite(const Gmpq &);
std::pair<double, double> to_interval (const Gmpq &);
#endif

#ifdef CGAL_USE_GMPXX
// GMPXX

CGAL_END_NAMESPACE
template <typename, typename> class __gmp_expr;
class __gmpz_value;
class __gmpq_value;
typedef __gmp_expr<__gmpz_value, __gmpz_value> mpz_class;
typedef __gmp_expr<__gmpq_value, __gmpq_value> mpq_class;
CGAL_BEGIN_NAMESPACE

template < typename T, typename U >
::__gmp_expr<T, T> sqrt(const ::__gmp_expr<T, U> &);

template < typename T, typename U >
double to_double(const ::__gmp_expr<T, U> &);

template < typename T, typename U >
bool is_finite(const ::__gmp_expr<T, U> &);

template < typename T, typename U >
bool is_valid(const ::__gmp_expr<T, U> &);

template < typename T, typename U >
std::pair<double,double> to_interval (const ::__gmp_expr<T, U> &);

std::pair<double, double> to_interval (const mpz_class &);

std::pair<double, double> to_interval (const mpq_class &);

template < typename T, typename U >
::__gmp_expr<T, T> abs(const ::__gmp_expr<T, U>&);

template < typename T, typename U >
::__gmp_expr<T, T> square(const ::__gmp_expr<T, U>&);


template < typename T, typename U >
Sign sign(const ::__gmp_expr<T, U> &);

template < typename T, typename U1, typename U2 >
Comparison_result
compare(const ::__gmp_expr<T, U1> &, const ::__gmp_expr<T, U2> &);

template < typename T, typename U >
bool is_zero(const ::__gmp_expr<T, U> &);

template < typename T, typename U >
bool is_one(const ::__gmp_expr<T, U> &);

template < typename T, typename U >
bool is_positive(const ::__gmp_expr<T, U> &);

template < typename T, typename U >
bool is_negative(const ::__gmp_expr<T, U> &);
#endif // CGAL_USE_GMPXX

// CORE::Expr
#ifdef CGAL_USE_CORE
CGAL_END_NAMESPACE

namespace CORE {
  class Expr;
}

CGAL_BEGIN_NAMESPACE

double to_double(const CORE::Expr &);
CORE::Expr sqrt(const CORE::Expr &);
bool is_finite(const CORE::Expr &);
bool is_valid(const CORE::Expr &);
Sign sign(const CORE::Expr&);
Comparison_result compare(const CORE::Expr&, const CORE::Expr&);
std::pair<double,double> to_interval (const CORE::Expr &);

CGAL_END_NAMESPACE
#endif

// specializations for Quotient

#include <CGAL/Quotient.h>

#ifdef CGAL_USE_GMP
#  include <CGAL/Gmpz.h>
#endif

#include <CGAL/MP_Float.h>

CGAL_BEGIN_NAMESPACE

#ifdef CGAL_USE_GMP
double to_double(const Quotient<Gmpz>&);
#endif

double to_double(const Quotient<MP_Float>&);
std::pair<double,double> to_interval(const Quotient<MP_Float>&);

CGAL_END_NAMESPACE

#endif // CGAL_CFG_HAS_TWO_STAGE_NAME_LOOKUP

#include <CGAL/number_utils_classes.h>

#endif // CGAL_NUMBER_TYPE_BASIC_H
