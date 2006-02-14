// Copyright (c) 1999-2005  Utrecht University (The Netherlands),
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

#ifndef CGAL_LAZY_EXACT_NT_FWD_H
#define CGAL_LAZY_EXACT_NT_FWD_H

// Forward declarations of functions on Lazy_exact_nt.

CGAL_BEGIN_NAMESPACE

template <typename ET> class Lazy_exact_nt;

template <typename ET>
double to_double(const Lazy_exact_nt<ET> &);

template <typename ET>
std::pair<double,double> to_interval(const Lazy_exact_nt<ET> &);

template <typename ET>
Sign sign(const Lazy_exact_nt<ET> &);

template <typename ET1, typename ET2>
Comparison_result
compare(const Lazy_exact_nt<ET1> &, const Lazy_exact_nt<ET2> &);

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

template <typename ET>
Lazy_exact_nt<ET> gcd(const Lazy_exact_nt<ET> &, const Lazy_exact_nt<ET> &);

CGAL_END_NAMESPACE

#endif // CGAL_LAZY_EXACT_NT_FWD_H
