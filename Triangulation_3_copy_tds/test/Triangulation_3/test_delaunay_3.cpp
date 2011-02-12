// Copyright (c) 1998  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
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
// Author(s)     : Francois Rebufat

#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Triangulation_cell_base_with_circumcenter_3.h>

bool del=true;

#include <CGAL/_test_types.h>
#include <CGAL/_test_cls_delaunay_3.h>

// Explicit instantiation of the whole class :
template class CGAL::Delaunay_triangulation_3<K>;

int main()
{
  typedef CGAL::Delaunay_triangulation_3<EPIC>  Cls;
  typedef CGAL::Delaunay_triangulation_3<EPEC>  Cls_with_epec;

  _test_cls_delaunay_3( Cls() );
  _test_cls_delaunay_3( Cls_with_epec() );

  // Second version for the circumcenter storing cell base class.
  typedef CGAL::Triangulation_vertex_base_3<K>                 Vb;
  typedef CGAL::Triangulation_cell_base_with_circumcenter_3<K> Cb;
  typedef CGAL::Triangulation_data_structure_3<Vb, Cb>         TDS;
  typedef CGAL::Delaunay_triangulation_3<K, TDS>               Cls_circumcenter;

  _test_cls_delaunay_3( Cls_circumcenter() );

  return 0;
}
