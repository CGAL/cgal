// Copyright (c) 1998  Utrecht University (The Netherlands),
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

// CGAL_CFG_NO_TMPL_IN_TMPL_DEPENDING_FUNCTION_PARAM.C
// ---------------------------------------------------------------------
// A short test program to evaluate a C++ compiler.
// This program is used by cgal_configure.
// The following documentation will be pasted in the generated configfile.
// ---------------------------------------------------------------------

//| G++ 2.95.2 has problems with member functions implemented outside of
//| the class body if this member function has a parameter type that is
//| dependant on a template in the template parameter list of the class. A
//| workaround would be to implement the member function inline in the class.
//| The following definition is set if this error error occurs.

template < template < class T> class T_HDS>
struct Container {
    typedef T_HDS<int>           HDS;
    typedef typename HDS::Handle Handle;
    void foo( Handle h);
    //void foo( Handle h) {} // workaround: implement foo inline.
};

template < template < class T> class T_HDS >
void Container<T_HDS>::foo( Handle ) {}

template <class T>
struct A {
    typedef int* Handle;
};

int main() {
    Container<A> c;
    c.foo(0);
    return 0;
}

// EOF //
