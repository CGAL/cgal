// Copyright (c) 1997-2001  Utrecht University (The Netherlands),
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

// CGAL_CFG_MATCHING_BUG_3.C
// ---------------------------------------------------------------------
// This program is used by cgal_configure.
// The following documentation will be pasted in the generated configfile.
// ---------------------------------------------------------------------

//| This flag is set, if the compiler does not match function arguments 
//| of pointer type correctly, when the return type depends on 
//| the parameter's type. (e.g. sun C++ 5.3)

template < class T > struct A     { typedef typename T::CCC CCC; };
template < class T > struct A<T*> { typedef typename T::CCC CCC; };

template < class T >
typename A< T >::CCC
foo(T)
{
  typedef typename A< T >::CCC C;
  return C();
}

struct B { typedef int CCC; };

int main()
{
  B *p = 0;
  return foo(p);
}

// EOF //
