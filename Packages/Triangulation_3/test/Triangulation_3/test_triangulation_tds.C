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
// source        :
// file          : test_triangulation_tds.C
// revision      : 
// revision_date : 
// author(s)     : Francois Rebufat
//
// coordinator   : INRIA Sophia-Antipolis
// ============================================================================

#include <CGAL/Triangulation_data_structure_3.h>

#include <CGAL/_test_cls_tds_3.C>

typedef CGAL::Triangulation_data_structure_3<>               Tds;

// Explicit instantiation :
// template class CGAL::Triangulation_data_structure_3<>;
// SunPRO requires putting all default arguments here.
template class CGAL::Triangulation_data_structure_3<
                     CGAL::Triangulation_ds_vertex_base_3<>,
                     CGAL::Triangulation_ds_cell_base_3<> >;


int main()
{
  _test_cls_tds_3(Tds());
  return 0;
}
