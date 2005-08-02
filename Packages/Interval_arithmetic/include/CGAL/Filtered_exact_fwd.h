// Copyright (c) 1998-2005  Utrecht University (The Netherlands),
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

#ifndef CGAL_FILTERED_EXACT_FWD_H
#define CGAL_FILTERED_EXACT_FWD_H

// Forward declarations

CGAL_BEGIN_NAMESPACE

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

CGAL_END_NAMESPACE

#endif // CGAL_FILTERED_EXACT_FWD_H
