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
// $URL: svn+ssh://fcacciola@scm.gforge.inria.fr/svn/cgal/trunk/Installation/config/support/test_BOOST.cpp $
// $Id: test_BOOST.cpp 32424 2006-07-12 09:26:22Z spion $
// 
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
