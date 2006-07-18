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

#ifndef CGAL_NEF_POLYNOMIAL_FWD_H
#define CGAL_NEF_POLYNOMIAL_FWD_H

#include <CGAL/Quotient_fwd.h>

// Forward declarations of functions over Polynomial and Nef_polynomial

CGAL_BEGIN_NAMESPACE

namespace Nef {
template <typename> class Polynomial;


template <typename ET>
double to_double(const Polynomial<ET> &);

//template <typename ET>
//std::pair<double,double> to_interval(const Polynomial<ET> &);

template <typename ET>
Sign sign(const Polynomial<ET> &);


template <typename ET>
Polynomial<ET> abs(const Polynomial<ET> &);

template <typename ET>
bool is_finite(const Polynomial<ET> &);

template <typename ET>
bool is_valid(const Polynomial<ET> &);

template <typename ET>
Polynomial<ET> gcd(const Polynomial<ET> &, const Polynomial<ET> &);

}
// Nef_polynomial

template <typename> class Nef_polynomial;

template <typename ET>
double to_double(const Nef_polynomial<ET> &);

template <class NT>
std::pair<double,double> to_interval(const Nef_polynomial<NT>& p);

template <typename ET>
Nef_polynomial<ET> gcd(const Nef_polynomial<ET> &, const Nef_polynomial<ET> &);

//template <typename ET>
//double to_double(const Quotient<Nef_polynomial<ET> >&);



using Nef::to_double;
using Nef::sign;
using Nef::abs;
using Nef::is_finite;
using Nef::is_valid;
using Nef::gcd;

CGAL_END_NAMESPACE

#endif // CGAL_NEF_POLYNOMIAL_FWD_H
