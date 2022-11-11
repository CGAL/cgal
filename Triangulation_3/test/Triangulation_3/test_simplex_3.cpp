// Copyright (c) 1998  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Nico Kruithof

//#define CGAL_TRIANGULATION_DONT_USE_SHORT_NAMES

#include <CGAL/Triangulation_3.h>

bool del = false;

#include <CGAL/_test_types.h>
#include <CGAL/_test_cls_triangulation_simplex_3.h>

// Explicit instantiation of the whole class :
template class CGAL::Triangulation_3<K>;

int main()
{
  typedef CGAL::Triangulation_3<K>                               Cls3;

  _test_cls_triangulation_simplex_3( Cls3() );

  return 0;
}
