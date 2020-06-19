// Copyright (c) 2005
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
