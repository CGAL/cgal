// Copyright (c) 2005  
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
// Author(s)     : Sylvain Pion

// Tests if BOOST is available.

#include <iostream>
#include <boost/version.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>

using boost::tuple;
using boost::make_tuple;
using boost::tie;
using boost::get;

int main()
{
    tuple<int, double> a, b, c;
    a = tuple<int, double>();
    b = tuple<int, double>(1);
    c = tuple<int, double>(1, 3.14);
    a = make_tuple(1, 2.57);

    int i;
    double d;
    tie(i, d) = a;
    i = a.get<0>();
    d = a.get<1>();

    std::cout << "version=" << BOOST_VERSION/100000 << "."
                            << ((BOOST_VERSION / 100) % 100) << "."
                            << BOOST_VERSION % 100 << std::endl;

    return 0;
}
