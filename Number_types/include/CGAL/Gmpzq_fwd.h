// Copyright (c) 1999,2003,2004,2005  Utrecht University (The Netherlands),
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
// $URL$
// $Id$
// 
//
// Author(s)     : Andreas Fabri, Stefan Schirra, Sylvain Pion

#ifndef CGAL_GMPZQ_FWD_H
#define CGAL_GMPZQ_FWD_H

#include <CGAL/Quotient_fwd.h>
#include <CGAL/Root_of_2_fwd.h>

// Forward declarations of functions over Gmpz and Gmpq.

CGAL_BEGIN_NAMESPACE

class Gmpz;

double to_double(const Gmpz&);
Sign sign(const Gmpz &);
bool is_valid(const Gmpz &);
bool is_finite(const Gmpz &);
Gmpz sqrt(const Gmpz &);
Gmpz div(const Gmpz &, const Gmpz &);
Gmpz gcd(const Gmpz &, const Gmpz &);
Gmpz gcd(const Gmpz &, int);
std::pair<double, double> to_interval (const Gmpz &);

class Gmpq;

double to_double(const Gmpq&);
Sign sign(const Gmpq &);
bool is_valid(const Gmpq &);
bool is_finite(const Gmpq &);
std::pair<double, double> to_interval (const Gmpq &);

double to_double(const Quotient<Gmpz>&);

/* FIX THE make_root_of_2 first
Root_of_2< Gmpz >
make_root_of_2(const Gmpq &a, const Gmpq &b,
               const Gmpq &c, bool d);

Root_of_2< Gmpz >
make_root_of_2(const Gmpz &a, const Gmpz &b,
               const Gmpz &c, bool d);
*/

class Gmpzf;

double to_double(const Gmpzf&);
Sign sign(const Gmpzf &);
bool is_valid(const Gmpzf &);
bool is_finite(const Gmpzf &);
Gmpzf sqrt(const Gmpzf &);
Gmpzf div(const Gmpzf &, const Gmpzf &);
Gmpzf gcd(const Gmpzf &, const Gmpzf &);
Gmpzf gcd(const Gmpzf &, int);
// ToDo:
// std::pair<double, double> to_interval (const Gmpzf &);
Gmpzf exact_division(const Gmpzf&, const Gmpzf &);
double to_double(const Quotient<Gmpzf>&);
Comparison_result compare (const Gmpzf&, const Gmpzf &);
CGAL_END_NAMESPACE

#endif // CGAL_GMPZQ_FWD_H
