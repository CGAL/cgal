// Copyright (c) 2004  
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved. 
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
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

#ifdef CGAL_USE_GMP
#  include <CGAL/Gmpz.h>
#  include <CGAL/Gmpq.h>
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
// #  include <CGAL/CORE_Expr.h>
namespace CORE {
class Expr;
}
#endif

namespace CGAL { namespace internal {

// A class which tells the prefered "exact number type" corresponding to a type.

// The default template chooses Gmpq or Quotient<MP_Float>.
// It should support the built-in types.
template < typename >
struct Exact_type_selector
#ifdef CGAL_USE_GMP
{ typedef Gmpq Type; };
#else
{ typedef Quotient<MP_Float> Type; };
#endif

template <>
struct Exact_type_selector<MP_Float>
{ typedef Quotient<MP_Float> Type; };

template <>
struct Exact_type_selector<Quotient<MP_Float> >
{ typedef Quotient<MP_Float> Type; };

// And we specialize for the following types :
#ifdef CGAL_USE_GMP
template <>
struct Exact_type_selector<Gmpz>
{ typedef Gmpq  Type; };

template <>
struct Exact_type_selector<Gmpq>
{ typedef Gmpq  Type; };
#endif

#ifdef CGAL_USE_GMPXX
template <>
struct Exact_type_selector< ::mpz_class>
{ typedef ::mpq_class  Type; };

template <>
struct Exact_type_selector< ::mpq_class>
{ typedef ::mpq_class  Type; };
#endif

#ifdef CGAL_USE_LEDA
template <>
struct Exact_type_selector<leda_integer>
{ typedef leda_rational  Type; };

template <>
struct Exact_type_selector<leda_rational>
{ typedef leda_rational  Type; };

template <>
struct Exact_type_selector<leda_real>
{ typedef leda_real  Type; };
#endif

#ifdef CGAL_USE_CORE
template <>
struct Exact_type_selector<CORE::Expr>
{ typedef CORE::Expr  Type; };
#endif

template < typename ET >
struct Exact_type_selector<Lazy_exact_nt<ET> >
{
  // We have a choice here :
  // - using ET gets rid of the DAG computation as well as redoing the interval
  // - using Lazy_exact_nt<ET> might use sharper intervals.
  typedef ET  Type;
  // typedef Lazy_exact_nt<ET>  Type;
};

} } // namespace CGAL::internal

#endif // CGAL_INTERNAL_EXACT_TYPE_SELECTOR_H
