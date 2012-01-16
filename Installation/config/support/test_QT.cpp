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
// Author(s)     : various

// Tests if QT is available.

#include <iostream>
#include <qglobal.h>
#include <qarray.h>

QArray<int> fib( int num )                  // returns fibonacci array
{
    ASSERT( num > 2 );
    QArray<int> f( num );                   // array of ints

    f[0] = f[1] = 1;                        // initialize first two numbers
    for ( int i=2; i<num; i++ )
        f[i] = f[i-1] + f[i-2];     

    return f;
}

int main()
{
    std::cout << "version=" << QT_VERSION_STR << std::endl;

    QArray<int> a = fib( 6 );               // get 6 first fibonaccis

    if(4==a.find(5))
	return 0;
    return 1;
}
