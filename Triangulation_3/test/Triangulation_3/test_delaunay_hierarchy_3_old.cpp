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
#include <CGAL/Triangulation_hierarchy_3.h>

bool del=true;

#include <CGAL/_test_types.h>
#include <CGAL/_test_cls_delaunay_3.h>

int main()
{
  typedef CGAL::Triangulation_vertex_base_3<K>             Vbb;
  typedef CGAL::Triangulation_hierarchy_vertex_base_3<Vbb> Vb;
  typedef CGAL::Delaunay_triangulation_cell_base_3<K>      Cb;
  typedef CGAL::Triangulation_data_structure_3<Vb, Cb>     Tds;
  typedef CGAL::Delaunay_triangulation_3<K,Tds>            Dt;
  typedef CGAL::Triangulation_hierarchy_3<Dt>              Dh;

  _test_cls_delaunay_3( Dh() );

  return 0;
}
