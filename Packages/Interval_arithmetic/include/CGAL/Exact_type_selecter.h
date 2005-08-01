// Copyright (c) 2004  Utrecht University (The Netherlands),
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
// Author(s)     : Sylvain Pion

#ifndef CGAL_EXACT_TYPE_SELECTER_H
#define CGAL_EXACT_TYPE_SELECTER_H

// This is an undocumented private helper for Filtered_kernel.

#include <CGAL/basic.h>
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
#  include <CGAL/CORE_Expr.h>
#endif

CGAL_BEGIN_NAMESPACE

// A class which tells the prefered exact number type corresponding to a type.

// The default template chooses Quotient<MP_Float>.
// It should support the built-in types, MP_Float, Quotient<MP_Float>.
template < typename >
struct Exact_type_selecter
{ typedef Quotient<MP_Float> Type; };

// And we specialize for the following types :
#ifdef CGAL_USE_GMP
template <>
struct Exact_type_selecter<Gmpz>
{ typedef Gmpq  Type; };

template <>
struct Exact_type_selecter<Gmpq>
{ typedef Gmpq  Type; };
#endif

#ifdef CGAL_USE_GMPXX
template <>
struct Exact_type_selecter< ::mpz_class>
{ typedef ::mpq_class  Type; };

template <>
struct Exact_type_selecter< ::mpq_class>
{ typedef ::mpq_class  Type; };
#endif

#ifdef CGAL_USE_LEDA
template <>
struct Exact_type_selecter<leda_integer>
{ typedef leda_rational  Type; };

template <>
struct Exact_type_selecter<leda_rational>
{ typedef leda_rational  Type; };

template <>
struct Exact_type_selecter<leda_real>
{ typedef leda_real  Type; };
#endif

#ifdef CGAL_USE_CORE
template <>
struct Exact_type_selecter<CORE::Expr>
{ typedef CORE::Expr  Type; };
#endif

template < typename ET >
struct Exact_type_selecter<Lazy_exact_nt<ET> >
{
  // We have a choice here :
  // - using ET gets rid of the DAG computation as well as redoing the interval
  // - using Lazy_exact_nt<ET> might use sharper intervals.
  typedef ET  Type;
  // typedef Lazy_exact_nt<ET>  Type;
};

CGAL_END_NAMESPACE

#endif // CGAL_EXACT_TYPE_SELECTER_H
