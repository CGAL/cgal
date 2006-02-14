// Copyright (c) 2004  Utrecht University (The Netherlands),
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
// Author(s)     : Menelaos Karavelas, Sylvain Pion

// ---------------------------------------------------------------------
// A short test program to evaluate a C++ compiler.
// This program is used by cgal_configure.
// The following documentation will be pasted in the generated configfile.
// ---------------------------------------------------------------------

//| This flag is set if the compiler doesn't support "using Base::Member;" in
//| a derived class.  The workaround is to write a forwarder.
//| At least g++ 2.95 and SunPro CC 5.3 have this bug.

template < typename T >
void kill_unused_warning(const T&) {}

struct A
{
  int operator()(int a) const { return a; }
};

template<class T>
struct B : public A
{
  using A::operator();

  T operator()(const T& t) const
  { return t+1; }
};


int main()
{
  B<double> b;
  double d1(0), d2(0);

  d2 = b(d1);

  kill_unused_warning(d1);
  kill_unused_warning(d2);

  return 0;
}
