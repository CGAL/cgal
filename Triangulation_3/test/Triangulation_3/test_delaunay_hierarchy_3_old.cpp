// Copyright (c) 1998,2001  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// 
//
// Author(s)     : Mariette Yvinec, Sylvain Pion

#define CGAL_NO_DEPRECATION_WARNINGS

#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Triangulation_hierarchy_3.h>

bool del=true;

#include <CGAL/_test_types.h>
#include <CGAL/_test_cls_delaunay_3.h>

int main()
{
  typedef CGAL::Triangulation_vertex_base_3<K>             Vbb;
  typedef CGAL::Triangulation_hierarchy_vertex_base_3<Vbb> Vb;
  typedef CGAL::Triangulation_data_structure_3<Vb>         Tds;
  typedef CGAL::Delaunay_triangulation_3<K,Tds>            Dt;
  typedef CGAL::Triangulation_hierarchy_3<Dt>              Dh;

  _test_cls_delaunay_3( Dh() );

  return 0;
}
