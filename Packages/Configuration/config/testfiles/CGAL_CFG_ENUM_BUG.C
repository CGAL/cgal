// Copyright (c) 1997  Utrecht University (The Netherlands),
// ETH Zurich (Switzerland), Freie Universitaet Berlin (Germany),
// INRIA Sophia-Antipolis (France), Martin-Luther-University Halle-Wittenberg
// (Germany), Max-Planck-Institute Saarbrucken (Germany), RISC Linz (Austria),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : various

// CGAL_CFG_ENUM_BUG.C
// ---------------------------------------------------------------------
// This program is used by cgal_configure.
// The following documentation will be pasted in the generated configfile.
// ---------------------------------------------------------------------

//| This flag is set, if the compiler does not promote enumeration types
//| (which depend on a template parameter) correctly when they are used 
//| as int template arguments. (e.g. Borland 5.5)

struct F {
  enum { a = 1 };
};

template < int i > struct B;
template <> struct B< 1 > {};

template < class T >
struct C {
  enum { ar = T::a };
  B< ar > b;
  // gives
  // Borland C++ 5.5.1 for Win32 Copyright (c) 1993, 2000 Borland
  // Error E2450 Undefined structure 'B<0>' in function main()
  
  // using 
  // B< T::a > b;
  // instead gives
  // Borland C++ 5.5.1 for Win32 Copyright (c) 1993, 2000 Borland
  // Error E2401 Invalid template argument list
  // Error E2040 Declaration terminated incorrectly
};

int main()
{
  C< F > c;
  (void) c;
  return 0;
}

// EOF //
