// Copyright (c) 1998  Utrecht University (The Netherlands),
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
// Author(s)     : various

// CGAL_CFG_NO_PARTIAL_CLASS_TEMPLATE_SPECIALISATION.C
// ---------------------------------------------------------------------

//| If a compiler doesn't support partial specialisation of class templates,
//| the flag CGAL_CFG_NO_PARTIAL_CLASS_TEMPLATE_SPECIALISATION is set.

template < class T >
struct B {
  T a;
};

template < class T >
struct C {
  T b;
};

template < class T >
struct C< B < T > > {
  int c;
};

int
main()
{
  C< int > t1;
  C< B< int > > t2;
  t1.b = 1;
  t2.c = 1;
  return 0;
}

// EOF //
