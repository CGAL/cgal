// Copyright (c) 2004
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Lutz Kettner, Sylvain Pion

//| A basic test for the STL.
//| If it fails, it probably means a bad CGAL installation.

#undef NDEBUG
#include <cassert>
#include <algorithm>
#include <list>
#include <iterator>

using std::list;

list<char> lst( const char* s)
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
