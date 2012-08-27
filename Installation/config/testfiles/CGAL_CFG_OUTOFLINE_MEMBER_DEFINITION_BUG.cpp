// Copyright (c) 2004  
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved. 
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
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
