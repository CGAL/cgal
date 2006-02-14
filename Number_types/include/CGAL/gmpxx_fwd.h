// Copyright (c) 2002,2003,2005  Utrecht University (The Netherlands),
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
 
#ifndef CGAL_GMPXX_FWD_H
#define CGAL_GMPXX_FWD_H

#ifdef CGAL_USE_GMPXX

// Forward declarations for the GMPXX types.

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

CGAL_END_NAMESPACE

#endif // CGAL_USE_GMPXX

#endif // CGAL_GMPXX_FWD_H
