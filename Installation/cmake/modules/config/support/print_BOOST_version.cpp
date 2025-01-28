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
// Author(s)     : Sylvain Pion

// Tests if BOOST is available.

#include <iostream>

#include <boost/version.hpp>

int main()
{
    std::cout << "version=" << BOOST_VERSION/100000 << "."
                            << ((BOOST_VERSION / 100) % 100) << "."
                            << BOOST_VERSION % 100 << std::endl;

    return 0;
}
