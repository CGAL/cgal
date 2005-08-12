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

#ifndef CGAL_SAME_UNCERTAINTY_H
#define CGAL_SAME_UNCERTAINTY_H

#include <CGAL/basic.h>

CGAL_BEGIN_NAMESPACE

template < typename T1, typename T2 >
struct Same_uncertainty
{
  typedef T1 type;
};

template < typename T > struct Uncertain;
template < typename T > struct Sgn;

template < typename T1, typename T2 >
struct Same_uncertainty < T1, Uncertain<T2> >
{
  typedef Uncertain<T1> type;
};

// Short cut to extract uncertainty from a number type directly.
template < typename T1, typename NT >
struct Same_uncertainty_nt
  : Same_uncertainty <T1, typename Sgn<NT>::result_type > {};

CGAL_END_NAMESPACE

#endif // CGAL_SAME_UNCERTAINTY_H
