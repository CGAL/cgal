// Copyright (c) 2007  INRIA Sophia-Antipolis (France).
// All rights reserved.
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
// $URL: svn+ssh://scm.gforge.inria.fr/svn/cgal/trunk/Installation/config/support/test_BOOST_PROGRAM_OPTIONS.cpp $
// $Id: test_BOOST_PROGRAM_OPTIONS.cpp 37323 2007-03-20 19:14:02Z reichel $
// 
//
// Author(s)     : Sylvain Pion

// Tests if BOOST.THREAD is available.

#include <iostream>
#include <boost/version.hpp>
#include <boost/thread/tss.hpp>

int main()
{
  std::cout << "version=" << BOOST_VERSION/100000 << "."
            << ((BOOST_VERSION / 100) % 100) << "."
            << BOOST_VERSION % 100 << std::endl;

  boost::thread_specific_ptr<int> z;
  if (z.get() == NULL) {
    z.reset(new int(1));
  }
  if(*z.get() == 1) {
    return 0;
  }
  else {
    return 1;
  }  
}
