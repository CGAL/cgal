// ============================================================================
//
// Copyright (c) 1998,2001 The CGAL Consortium
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
// file          : test/Triangulation3/test_delaunay_hierarchy_3.C
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Mariette Yvinec, Sylvain Pion
//
// coordinator   : INRIA Sophia-Antipolis
// ============================================================================

#include <CGAL/basic.h>
#include <iostream>

#include <CGAL/_test_types.h>
#include <CGAL/Simple_cartesian.h>

#include <CGAL/Point_3.h>
#include <CGAL/Segment_3.h>
#include <CGAL/Triangle_3.h>
#include <CGAL/Line_3.h>
#include <CGAL/Direction_3.h>
#include <CGAL/Ray_3.h>

#include <CGAL/Triangulation_geom_traits_3.h>
#include <CGAL/Triangulation_data_structure_3.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Triangulation_hierarchy_3.h>

#include <CGAL/_test_cls_delaunay_3.C>

bool del=true;

int main()
{
  std::cout << "Testing Delaunay_hierarchy_3 " << std::endl; 
  std::cout << " with Cartesian<double> points "<<  std::endl;

  typedef double                                           Coord_type;
  typedef CGAL::Simple_cartesian<Coord_type>               Rep;
  // typedef CGAL::Triangulation_geom_traits_3<Rep>          Gt;
  typedef Rep                                              Gt;
  typedef CGAL::Triangulation_vertex_base_3<Gt>            Vbb;
  typedef CGAL::Triangulation_hierarchy_vertex_base_3<Vbb> Vb;
  typedef CGAL::Triangulation_cell_base_3<Gt>              Ce;
  typedef CGAL::Triangulation_data_structure_3<Vb,Ce>      Tdsul;
  typedef CGAL::Delaunay_triangulation_3<Gt,Tdsul>         Dt;
  typedef CGAL::Triangulation_hierarchy_3<Dt>              Dh;

  _test_cls_delaunay_3( Dh() );

  return 0;
}
