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
#include <CGAL/float.h>
#include <CGAL/double.h>
#include <CGAL/long_double.h>
#include <CGAL/int.h>

// Including all number type files is necessary for compilers implementing
// two-stage name lookup (like g++ >= 3.4).
// A nicer solution needs more thought.

#ifdef CGAL_CFG_HAS_TWO_STAGE_NAME_LOOKUP

#include <CGAL/Quotient_fwd.h>
#include <CGAL/Interval_nt_fwd.h>
#include <CGAL/Lazy_exact_nt_fwd.h>
#include <CGAL/Gmpzq_fwd.h>
#include <CGAL/gmpxx_fwd.h>
#include <CGAL/Number_type_checker_fwd.h>

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


// Polynomial

template <typename> class Polynomial;

template <typename ET>
double to_double(const Polynomial<ET> &);

template <typename ET>
std::pair<double,double> to_interval(const Polynomial<ET> &);

template <typename ET>
Sign sign(const Polynomial<ET> &);


template <typename ET>
Polynomial<ET> abs(const Polynomial<ET> &);

template <typename ET>
bool is_finite(const Polynomial<ET> &);

template <typename ET>
bool is_valid(const Polynomial<ET> &);

template <typename ET>
Polynomial<ET> gcd(const Polynomial<ET> &, const Polynomial<ET> &);

// Nef_polynomial

template <typename> class Nef_polynomial;

template <typename ET>
double to_double(const Nef_polynomial<ET> &);

template <typename ET>
Nef_polynomial<ET> gcd(const Nef_polynomial<ET> &, const Nef_polynomial<ET> &);


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
#endif // CGAL_USE_CORE


// specializations for Quotient

CGAL_END_NAMESPACE
#include <CGAL/Quotient.h>

#ifdef CGAL_USE_GMP
#  include <CGAL/Gmpz.h>
#endif // CGAL_USE_GMP

#include <CGAL/MP_Float.h>

CGAL_BEGIN_NAMESPACE

#ifdef CGAL_USE_GMP
double to_double(const Quotient<Gmpz>&);
#endif // CGAL_USE_GMP

double to_double(const Quotient<MP_Float>&);
std::pair<double,double> to_interval(const Quotient<MP_Float>&);

CGAL_END_NAMESPACE

#endif // CGAL_CFG_HAS_TWO_STAGE_NAME_LOOKUP

#include <CGAL/number_utils_classes.h>

#endif // CGAL_NUMBER_TYPE_BASIC_H
