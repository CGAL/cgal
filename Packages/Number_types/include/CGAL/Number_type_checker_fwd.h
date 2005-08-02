// Copyright (c) 2005  Utrecht University (The Netherlands),
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

#ifndef CGAL_NUMBER_TYPE_CHECKER_FWD_H
#define CGAL_NUMBER_TYPE_CHECKER_FWD_H

// Forward declarations

CGAL_BEGIN_NAMESPACE

template < typename NT1, typename NT2, typename Cmp >
class Number_type_checker;

template < typename NT1, typename NT2, typename Cmp >
double
to_double(const Number_type_checker<NT1, NT2, Cmp> &);

template < typename NT1, typename NT2, typename Cmp >
std::pair<double, double>
to_interval(const Number_type_checker<NT1, NT2, Cmp> &);

template < typename NT1, typename NT2, typename Cmp >
Number_type_checker<NT1, NT2, Cmp>
sqrt(const Number_type_checker<NT1, NT2, Cmp> &);

template < typename NT1, typename NT2, typename Cmp >
bool
is_finite(const Number_type_checker<NT1, NT2, Cmp> &);

template < typename NT1, typename NT2, typename Cmp >
bool
is_valid(const Number_type_checker<NT1, NT2, Cmp> &);

template < typename NT1, typename NT2, typename Cmp >
Sign
sign(const Number_type_checker<NT1, NT2, Cmp> &);

template < typename NT1, typename NT2, typename Cmp >
Comparison_result
compare(const Number_type_checker<NT1, NT2, Cmp> &,
        const Number_type_checker<NT1, NT2, Cmp> &);

CGAL_END_NAMESPACE

#endif // CGAL_NUMBER_TYPE_CHECKER_FWD_H
