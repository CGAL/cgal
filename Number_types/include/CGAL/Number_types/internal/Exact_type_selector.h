// Copyright (c) 2004
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Sylvain Pion

#ifndef CGAL_INTERNAL_EXACT_TYPE_SELECTOR_H
#define CGAL_INTERNAL_EXACT_TYPE_SELECTOR_H

// This is an undocumented private helper for Filtered_kernel.

#include <CGAL/number_type_basic.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Quotient.h>
#include <CGAL/Lazy_exact_nt.h>

#include <CGAL/boost_mp.h>

#ifdef CGAL_USE_GMP
#  include <CGAL/Gmpz.h>
#  include <CGAL/Gmpq.h>
#  include <CGAL/Gmpzf.h>
#  include <CGAL/Mpzf.h>
#endif
#ifdef CGAL_USE_GMPXX
#  include <CGAL/gmpxx.h>
#endif
#ifdef CGAL_USE_LEDA
#  include <CGAL/leda_integer.h>
#  include <CGAL/leda_rational.h>
#  include <CGAL/leda_real.h>
#endif
#ifdef CGAL_USE_CORE
// # include <CGAL/CORE_Expr.h>
namespace CORE {
class Expr;
}
#endif

namespace CGAL { namespace internal {

// Two classes which tell the preferred "exact number types" corresponding to a type.
// Exact_ring_selector<double> and Exact_field_selector<double> are used by EPICK as exact number type
// to answer predicates at the end of the filtering chain of predicates and EPECK uses
// Exact_field_selector<double> for as its exact number type.


enum ENT_backend_choice
{
  GMP_BACKEND,
  GMPXX_BACKEND,
  BOOST_GMP_BACKEND,
  BOOST_BACKEND,
  LEDA_BACKEND,
  MP_FLOAT_BACKEND
};

template <ENT_backend_choice>
struct Exact_NT_backend;

#ifdef CGAL_USE_GMP
template <>
struct Exact_NT_backend<GMP_BACKEND>
{
  typedef Gmpq Rational;
#ifdef CGAL_HAS_MPZF
  typedef Mpzf Integer;
#else
  typedef Gmpzf Integer;
#endif
};
#endif

#ifdef CGAL_USE_GMPXX
template <>
struct Exact_NT_backend<GMPXX_BACKEND>
{
  typedef mpq_class Rational;
  typedef Exact_NT_backend<GMP_BACKEND>::Integer Integer;
};
#endif

#ifdef CGAL_USE_BOOST_MP
template <>
struct Exact_NT_backend<BOOST_GMP_BACKEND>
{
  typedef BOOST_gmp_arithmetic_kernel::Rational Rational;
  typedef Exact_NT_backend<GMP_BACKEND>::Integer Integer;
};

template <>
struct Exact_NT_backend<BOOST_BACKEND>
{
// See the discussion in https://github.com/CGAL/cgal/pull/3614
// This is disabled for now because cpp_rational is even slower than Quotient<MP_Float>. Quotient<cpp_int> will be a good candidate after some polishing.
// In fact, the new version of cpp_rational from here: https://github.com/boostorg/multiprecision/pull/366
// is much better than Quotient<cpp_int> because it is using smart gcd and is well-supported
// while Quotient does not. Though, we can still use it if needed.
#if BOOST_VERSION <= 107800
// See this comment: https://github.com/CGAL/cgal/pull/5937#discussion_r721533675
  typedef Quotient<boost::multiprecision::cpp_int> Rational;
#else
  typedef BOOST_cpp_arithmetic_kernel::Rational Rational;
#endif
  typedef cpp_float Integer;
};
#endif

#ifdef CGAL_USE_LEDA
template <>
struct Exact_NT_backend<LEDA_BACKEND>
{
  typedef leda_rational Rational;
  typedef leda_integer Integer;
};
#endif

template <>
struct Exact_NT_backend<MP_FLOAT_BACKEND>
{
  typedef Quotient<MP_Float> Rational;
  typedef MP_Float Integer;
};

constexpr ENT_backend_choice Default_Exact_nt_back_end =
#if BOOST_VERSION > 107900 && defined(CGAL_USE_BOOST_MP)
  BOOST_BACKEND;
#else // BOOST_VERSION > 107900
  #ifdef CGAL_USE_GMPXX
    GMPXX_BACKEND;
  #elif defined(CGAL_USE_GMP)
    #if defined(CGAL_USE_BOOST_MP)
      BOOST_GMP_BACKEND;
    #else
      GMP_BACKEND;
    #endif
  #elif defined(CGAL_USE_LEDA)
    LEDA_BACKEND;
  #else
    MP_FLOAT_BACKEND;
  #endif
#endif // BOOST_VERSION > 107900

template < typename >
struct Exact_field_selector
{
  using Type = typename Exact_NT_backend<Default_Exact_nt_back_end>::Rational;
};

// By default, a field is a safe choice of ring.
template < typename T >
struct Exact_ring_selector : Exact_field_selector < T > { };

template <>
struct Exact_ring_selector<double>
{
  using Type = typename Exact_NT_backend<Default_Exact_nt_back_end>::Integer;
};

template <>
struct Exact_ring_selector<float> : Exact_ring_selector<double> { };

template <>
struct Exact_field_selector<MP_Float>
{ typedef Quotient<MP_Float> Type; };

template <>
struct Exact_ring_selector<MP_Float>
{ typedef MP_Float Type; };

template <>
struct Exact_field_selector<Quotient<MP_Float> >
{ typedef Quotient<MP_Float> Type; };

// And we specialize for the following types :
#ifdef CGAL_USE_GMP
template <>
struct Exact_field_selector<Gmpz>
{ typedef Gmpq  Type; };

template <>
struct Exact_ring_selector<Gmpz>
{ typedef Gmpz  Type; };

template <>
struct Exact_ring_selector<Gmpzf>
{ typedef Gmpzf Type; };

template <>
struct Exact_field_selector<Gmpq>
{ typedef Gmpq  Type; };
#endif

#ifdef CGAL_USE_GMPXX
template <>
struct Exact_field_selector< ::mpz_class>
{ typedef ::mpq_class  Type; };

template <>
struct Exact_ring_selector< ::mpz_class>
{ typedef ::mpz_class  Type; };

template <>
struct Exact_field_selector< ::mpq_class>
{ typedef ::mpq_class  Type; };
#endif

#ifdef CGAL_USE_LEDA
template <>
struct Exact_field_selector<leda_integer>
{ typedef leda_rational  Type; };

template <>
struct Exact_ring_selector<leda_integer>
{ typedef leda_integer   Type; };

template <>
struct Exact_field_selector<leda_rational>
{ typedef leda_rational  Type; };

template <>
struct Exact_field_selector<leda_real>
{ typedef leda_real  Type; };
#endif

#ifdef CGAL_USE_CORE
template <>
struct Exact_field_selector<CORE::Expr>
{ typedef CORE::Expr  Type; };
#endif

template < typename ET >
struct Exact_field_selector<Lazy_exact_nt<ET> >
: Exact_field_selector<ET>
{
  // We have a choice here :
  // - using ET gets rid of the DAG computation as well as redoing the interval
  // - using Lazy_exact_nt<ET> might use sharper intervals.
  // typedef ET  Type;
  // typedef Lazy_exact_nt<ET>  Type;
};
template < typename ET >
struct Exact_ring_selector<Lazy_exact_nt<ET> >
: Exact_ring_selector<ET>
{};

#ifndef CGAL_NO_DEPRECATED_CODE
// Added for backward compatibility
template < typename ET >
struct Exact_type_selector : Exact_field_selector< ET > {};
#endif

} } // namespace CGAL::internal

#endif // CGAL_INTERNAL_EXACT_TYPE_SELECTOR_H
