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
  typedef CGAL::Triangulation_data_structure_3<CGAL::VertexWithPoint<EPIC>, CGAL::Cell4Delaunay<EPIC>, CGAL::Sequential_tag, CGAL::Index_tag> Tds_index;
  typedef CGAL::Triangulation_data_structure_3<CGAL::VertexWithPoint<EPEC>, CGAL::Cell4Delaunay<EPEC>, CGAL::Sequential_tag, CGAL::Index_tag> Tds_with_epec_index;
  typedef CGAL::Delaunay_triangulation_3<EPIC, CGAL::Fast_location> Dh;
  typedef CGAL::Delaunay_triangulation_3<EPIC, Tds_index, CGAL::Fast_location> Dh_index;
  typedef CGAL::Delaunay_triangulation_3<EPEC, CGAL::Fast_location> Dh_with_epec;
  typedef CGAL::Delaunay_triangulation_3<EPEC, Tds_with_epec_index, CGAL::Fast_location> Dh_with_epec_index;

  _test_cls_delaunay_3( Dh() );
  _test_cls_delaunay_3( Dh_with_epec() );
  _test_cls_delaunay_3( Dh_index() );
  _test_cls_delaunay_3( Dh_with_epec_index() );

  return 0;
}
