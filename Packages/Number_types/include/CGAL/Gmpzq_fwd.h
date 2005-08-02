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
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Andreas Fabri, Stefan Schirra, Sylvain Pion

#ifndef CGAL_GMPZQ_FWD_H
#define CGAL_GMPZQ_FWD_H

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

CGAL_END_NAMESPACE

#endif // CGAL_GMPZQ_FWD_H
