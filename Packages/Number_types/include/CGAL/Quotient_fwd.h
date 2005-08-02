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
// Author(s)     : Stefan Schirra, Sylvain Pion

#ifndef CGAL_QUOTIENT_FWD_H
#define CGAL_QUOTIENT_FWD_H

CGAL_BEGIN_NAMESPACE

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

CGAL_END_NAMESPACE

#include <CGAL/MP_Float.h>
#include <CGAL/Gmpz.h>

CGAL_BEGIN_NAMESPACE

double to_double(const Quotient<MP_Float>&);
std::pair<double,double> to_interval(const Quotient<MP_Float>&);

double to_double(const Quotient<Gmpz>&);

CGAL_END_NAMESPACE

#endif  // CGAL_QUOTIENT_FWD_H
