// Copyright (c) 1997  Utrecht University (The Netherlands),
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

// CGAL_CFG_EARLY_INSTANTIATION_BUG.C
// ---------------------------------------------------------------------
// A short test program to evaluate a C++ compiler / STL implementation
// whether it supports the new standard headers (i.e. without the .h suffix)
// This program is used by cgal_configure.
// The following documentation will be pasted in the generated configfile.
// ---------------------------------------------------------------------

//| This flag is set, if a compiler instantiates a class template A<>
//| too early, where A<> is declared in the scope of another class
//| template B<>, such that B<> is the template parameter of A<>.
//| This bug shows up on VC7.0.

template < class Container > class Iterator;

template < class A >
struct Container {
  typedef Container<A>    Self;
  typedef Iterator<Self>  iterator;

  // VC++ 7.0 instantiates iterator too early, unless passed by reference.
  // Solution: void erase(iterator&){}
  void erase(iterator) {}
};

template < class Container >
class Iterator {
  typedef typename Container::Self pipo;
};

int main() {
  Container<int> C;
  typedef Container<int>::iterator It;
  It i;
  C.erase(i);
  return 0;
}
