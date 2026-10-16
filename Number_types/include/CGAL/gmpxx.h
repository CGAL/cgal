// Copyright (c) 2002,2003
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

// This file gathers the necessary adaptors so that the following
// C++ number types that come with GMP can be used by CGAL :
// - mpz_class (see #include <CGAL/mpz_class.h>)
// - mpq_class (see #include <CGAL/mpq_class.h>)

// - mpf_class support is commented out until to_interval() is implemented.
//   It is probably not very useful with CGAL anyway.

// Note that GMP++ use the expression template mechanism, which makes things
// a little bit complicated in order to make square(x+y) work for example.
// Reading gmpxx.h shows that ::__gmp_expr<T, T> is the mp[zqf]_class proper,
// while ::__gmp_expr<T, U> is the others "expressions".

#ifndef CGAL_GMPXX_H
#define CGAL_GMPXX_H

#include <cstring> // needed by GMP 4.1.4 since <gmpxx.h> misses it.
#include <gmpxx.h>
#include <utility>

#include <CGAL/mpz_class.h>
#include <CGAL/mpq_class.h>
#include <CGAL/gmpxx_coercion_traits.h>

namespace CGAL {


template < typename T, typename U >
class Algebraic_structure_traits< ::__gmp_expr<T,U> >
    : public Algebraic_structure_traits< ::__gmp_expr<T,T> >{};

template < typename T, typename U >
class Real_embeddable_traits< ::__gmp_expr<T,U> >
    : public Real_embeddable_traits< ::__gmp_expr<T,T> >{};


template <class NT1, class NT2> struct NT_converter;

namespace internal {

template <class Integer, class Target>
struct Gmpxx_long_long_converter
{
  typedef Integer argument_type;
  typedef Target result_type;

  Target operator()(Integer value) const
  {
    unsigned long long magnitude = static_cast<unsigned long long>(value);
    if(value < 0)
      magnitude = 0 - magnitude;
    mpz_class integer;
    mpz_import(integer.get_mpz_t(), 1, 1, sizeof(magnitude), 0, 0, &magnitude);
    if(value < 0)
      mpz_neg(integer.get_mpz_t(), integer.get_mpz_t());
    return Target(integer);
  }
};

} // namespace internal

template <>
struct NT_converter<long long, mpz_class>
  : internal::Gmpxx_long_long_converter<long long, mpz_class> {};

template <>
struct NT_converter<unsigned long long, mpz_class>
  : internal::Gmpxx_long_long_converter<unsigned long long, mpz_class> {};

template <>
struct NT_converter<long long, mpq_class>
  : internal::Gmpxx_long_long_converter<long long, mpq_class> {};

template <>
struct NT_converter<unsigned long long, mpq_class>
  : internal::Gmpxx_long_long_converter<unsigned long long, mpq_class> {};

} //namespace CGAL

#include <CGAL/GMPXX_arithmetic_kernel.h>

#endif // CGAL_GMPXX_H
