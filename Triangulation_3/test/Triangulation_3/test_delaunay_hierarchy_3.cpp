// Copyright (c) 1998,2001  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Mariette Yvinec, Sylvain Pion

#include <CGAL/Delaunay_triangulation_3.h>

bool del=true;

#include <CGAL/_test_types.h>
#include <CGAL/_test_cls_delaunay_3.h>

int main()
{
  typedef CGAL::Delaunay_triangulation_3<EPIC, CGAL::Fast_location> Dh;
  typedef CGAL::Delaunay_triangulation_3<EPEC, CGAL::Fast_location> Dh_with_epec;

  _test_cls_delaunay_3( Dh() );
  _test_cls_delaunay_3( Dh_with_epec() );

  return 0;
}
