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

//#include <CGAL/Fixed_precision_nt.h> 
//#define Fixed_precision_nt FIX
#include <CGAL/Cartesian.h>

#include <CGAL/Triangulation_short_names_2.h>
#include <CGAL/Delaunay_hierarchy_2.h>

#include <CGAL/_test_cls_triangulation_short_2.C>
#include <CGAL/_test_cls_delaunay_hierarchy_2.C>

int main()
{
  std::cout << "Testing Delaunay_hierarchy_2 " << std::endl; 
  std::cout << " with Cartesian points, Fixed coordinates: "<<
    std::endl;
  // typedef CGAL::Fixed_precision_nt    Coord_type;
  typedef double                      Coord_type;
  typedef CGAL::Cartesian<Coord_type> Gt;
  typedef CGAL::Delaunay_hierarchy_vertex_base_2<Gt> Vb;
  typedef CGAL::Triangulation_face_base_2<Gt>  Fb;
  typedef CGAL::Triangulation_default_data_structure_2<Gt,Vb,Fb> Tds;
  typedef CGAL::Delaunay_hierarchy_2<Gt,Tds>  Dh;

 _test_cls_delaunay_hierarchy_2( Dh() );
}
