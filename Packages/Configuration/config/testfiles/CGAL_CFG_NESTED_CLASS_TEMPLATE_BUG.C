// Copyright (c) 2003  Utrecht University (The Netherlands),
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
// Author(s)     : Michael Hoffmann <hoffmann@inf.ethz.ch>

// CGAL_CFG_NESTED_CLASS_TEMPLATE_BUG.C
// ---------------------------------------------------------------------
// A short test program to evaluate a C++ compiler.
// This program is used by cgal_configure.
// The following documentation will be pasted in the generated configfile.
// ---------------------------------------------------------------------

//| This flag ist set, if a compiler cannot parse an instantiation of
//| a nested class template B which is defined in another class
//| template A. This bug shows up on VC7. A workaround is to pass the
//| instantiation of A to a helper class which extracts B. (see
//| CGAL/Kernel_d/Iso_box_d.h)

struct II {};

template < bool b > struct A;
template <> struct A<true> { 
  template <typename T> struct B {
    typedef int type;
  };
};
template <> struct A<false> { 
  template <typename T> struct B {
    typedef II type;
  };
};

template <bool b>
struct C {
  typedef typename A<b>::template B<int>::type D;
};

int main()
{
  C<true>::D i = 0;
  ++i;
  return 0;
}
