// Copyright (c) 2002  Utrecht University (The Netherlands),
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
// Author(s)     : Sylvain Pion

// CGAL_CFG_VC7_PRIVATE_TYPE_BUG.C
// ---------------------------------------------------------------------
// A short test program to evaluate a C++ compiler.
// This program is used by cgal_configure.
// The following documentation will be pasted in the generated configfile.
// ---------------------------------------------------------------------

//| This is a test-case for a bug in VC++ 7.0 beta2 that occurs in the kernel.
//| When the bug is present, CGAL_CFG_VC7_PRIVATE_TYPE_BUG is set.

template <class R>
class A
{
  typedef typename R::B B;
public:
 
  A() :i(0) {}
 
  int i;
};
 
template <class R>
class B
: public R::A
{
public:
  typedef typename R::A base;
  B() : base() {}
};
 
template <class FT>
class RR
{
public:
  typedef ::A<RR> A;
  typedef ::B<RR> B;
};
 
int main()
{
  typedef RR<int> r;
  A<r> a;
  B<r> b;
  (void) a;
  (void) b;
  return 0;
}
