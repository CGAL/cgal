// ============================================================================
//
// Copyright (c) 1998 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
// 
// release       :
// release_date  :
// 
// file          : test/Triangulation/test_delaunay_hierarchy_2.C
// source        : $RCSfile$
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Mariette Yvinec (Mariette.Yvinec@sophia.inria.fr)
//
// coordinator   : INRIA Sophia-Antipolis
// ============================================================================
#include <CGAL/basic.h>
#include <iostream>

#include <CGAL/_test_types.h>
#include <CGAL/Cartesian.h>

#include <CGAL/Point_2.h>
#include <CGAL/Segment_2.h>
#include <CGAL/Triangle_2.h>
#include <CGAL/Line_2.h>
#include <CGAL/Direction_2.h>
#include <CGAL/Ray_2.h>

#include <CGAL/Triangulation_short_names_2.h>
#include <CGAL/Triangulation_data_structure_using_list_2.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_hierarchy_2.h>

#include <CGAL/_test_cls_delaunay_hierarchy_2.C>

int main()
{
  std::cout << "Testing Delaunay_hierarchy_2 " << std::endl; 
  std::cout << " with Cartesian<double> points "<<  std::endl;

  typedef double                      Coord_type;
  typedef CGAL::Cartesian<Coord_type> Gt;
  typedef CGAL::Triangulation_vertex_base_2<Gt> Vbb;
  typedef CGAL::Triangulation_hierarchy_vertex_base_2<Vbb> Vb;
  typedef CGAL::Triangulation_face_base_2<Gt>  Fb;
  typedef CGAL::Triangulation_data_structure_using_list_2<Vb,Fb>  Tdsul;
  typedef CGAL::Delaunay_triangulation_2<Gt,Tdsul>  Dt;
  typedef CGAL::Triangulation_hierarchy_2<Dt>  Dh;

  _test_cls_delaunay_hierarchy_2( Dh() );

   return 0;
}
