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
// Author(s)     : Laurent Rineau

// CGAL_CFG_OUTOFLINE_TEMPLATE_MEMBER_DEFINITION_BUG_2.C
// ---------------------------------------------------------------------
// A short test program to evaluate a C++ compiler.
// This program is used by cgal_configure.
// The following documentation will be pasted in the generated configfile.
// ---------------------------------------------------------------------

//| The flag CGAL_CFG_OUTOFLINE_TEMPLATE_MEMBER_DEFINITION_BUG_2 is set,
//| if a compiler does not support the definition of member templates
//| out of line, i.e. outside class scope. The solution is to put the
//| definition inside the class.
//| The difference with CGAL_CFG_OUTOFLINE_TEMPLATE_MEMBER_DEFINITION_BUG
//| is the return type of the member template.SunPro 5.5 should be OK with
//| this code.

template <class A>
class B {
public:
  template<class C>
    void fct(C *i);
};

template<class A>
template<class C>
void
B<A>::fct(C *i){
}

int main(){
  return 0;
}
