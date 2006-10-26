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
// $URL$
// $Id$
// 
//
// Author(s)     : Sylvain Pion
 
#ifndef CGAL_GMPXX_FWD_H
#define CGAL_GMPXX_FWD_H

#ifdef CGAL_USE_GMPXX

#include <gmpxx.h>
#include <CGAL/Root_of_2_fwd.h>

CGAL_BEGIN_NAMESPACE

/* FIX gmpxx.h first
   template < typename T, typename U1, typename U2, typename U3 >
   typename Root_of_traits< ::__gmp_expr<T, T> >::RootOf_2
   make_root_of_2(const ::__gmp_expr< T, U1> & a,
   const ::__gmp_expr< T, U2> & b,
   const ::__gmp_expr< T, U3> & c,
   bool d);*/

CGAL_END_NAMESPACE

#endif // CGAL_USE_GMPXX

#endif // CGAL_GMPXX_FWD_H
