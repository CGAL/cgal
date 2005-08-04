// Copyright (c) 1999,2001  Utrecht University (The Netherlands),
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
// Author(s)     : Andreas Fabri

#ifndef CGAL_SIXTUPLE_H
#define CGAL_SIXTUPLE_H

CGAL_BEGIN_NAMESPACE

template < class T >
struct Sixtuple
{
  typedef T value_type;

  T  e0, e1, e2, e3, e4, e5;

  Sixtuple()
  {}

  Sixtuple(const T & a0, const T & a1, const T & a2,
           const T & a3, const T & a4, const T & a5)
    : e0(a0), e1(a1), e2(a2), e3(a3), e4(a4), e5(a5)
  {}
};

CGAL_END_NAMESPACE

#endif // CGAL_SIXTUPLE_H
