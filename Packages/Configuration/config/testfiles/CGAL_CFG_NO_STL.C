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
// Author(s)     : Lutz Kettner, Sylvain Pion

// CGAL_CFG_NO_STL.C
// ---------------------------------------------------------------------
// A short test program to evaluate a machine architecture.
// This program is used by cgal_configure.
// The following documentation will be pasted in the generated configfile.
// ---------------------------------------------------------------------

//| A basic test for the STL.
//| If it fails, it probably means a bad CGAL installation.

#include <algorithm>
#include <list>
#include <iterator>
#include <cassert>

using std::list;

list<char> lst( char* s)
{
    list<char> x;
    while (*s != '\0') x.push_back( *s++);
    return x;
}

int main()
{
    list<char> list1 = lst( "mark twain");
    std::reverse( list1.begin(), list1.end());
    assert( list1 == lst( "niawt kram"));
    return 0;
}
