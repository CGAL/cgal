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
// Author(s)     : Marc Glisse

//| This flag is set if the compiler complains about an ambiguity between
//| a type and itself when some members are defined out of line. This is
//| a Sun CC bug.

template < typename U >
struct B
{
    struct A
    {
        typedef char C;
    };
    A* f ( typename A :: C );
};

template < class U >
    typename B < U > :: A*
B < U > :: f ( typename A :: C )
{
    return 0;
}

int main()
{
    return 0;
}
