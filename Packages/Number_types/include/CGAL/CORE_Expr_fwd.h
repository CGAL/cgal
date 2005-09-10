// Copyright (c) 2002-2005  Utrecht University (The Netherlands),
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
 
#ifndef CGAL_CORE_EXPR_FWD_H
#define CGAL_CORE_EXPR_FWD_H

#ifdef CGAL_USE_CORE

// Forward declarations

namespace CORE
{
  class Expr;
}

CGAL_BEGIN_NAMESPACE

double to_double(const CORE::Expr &);
CORE::Expr sqrt(const CORE::Expr &);
bool is_finite(const CORE::Expr &);
bool is_valid(const CORE::Expr &);
Sign sign(const CORE::Expr&);
Comparison_result compare(const CORE::Expr&, const CORE::Expr&);
std::pair<double,double> to_interval (const CORE::Expr &);
io_Operator io_tag(const CORE::Expr &);

CGAL_END_NAMESPACE

#endif // CGAL_USE_CORE

#endif // CGAL_CORE_EXPR_FWD_H
