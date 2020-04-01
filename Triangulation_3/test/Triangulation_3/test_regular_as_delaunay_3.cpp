// Copyright (c) 2004  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Sylvain Pion

#include <CGAL/Regular_triangulation_3.h>


bool del=true;

#include "include/CGAL/_test_types.h"
#include "include/CGAL/_test_cls_delaunay_3.h"


int main()
{
  typedef CGAL::Regular_triangulation_3<K>    Cls;

  _test_cls_delaunay_3( Cls() );

  return 0;
}

// MipsPro prefers this after the other instantiations...
// Explicit instantiation of the whole class :
template class CGAL::Regular_triangulation_3<K>;
