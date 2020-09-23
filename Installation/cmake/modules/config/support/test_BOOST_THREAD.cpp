// Copyright (c) 2007  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
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
