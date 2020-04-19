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


#include <iostream>
#include <LEDA/system/basic.h>


int main()
{
  std::cout << "version=" << (  __LEDA__  / 100 ) << "." << ( __LEDA__ % 100 ) << std::endl;
  return 0;
}
