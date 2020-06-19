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
// Author(s)     : Francois Rebufat

#include <CGAL/internal/disable_deprecation_warnings_and_errors.h>

#include "test_dependencies.h"

#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Triangulation_data_structure_3.h>
#include <CGAL/Triangulation_vertex_base_3.h>
#include <CGAL/Delaunay_triangulation_cell_base_3.h>
#include <CGAL/Delaunay_triangulation_cell_base_with_circumcenter_3.h>

bool del=true;

#include <CGAL/_test_types.h>
#include <CGAL/_test_cls_delaunay_3.h>
#include <CGAL/_test_cls_parallel_triangulation_3.h>

// Explicit instantiation of the whole class :
template class CGAL::Delaunay_triangulation_3<K>;

int main()
{
  typedef CGAL::Delaunay_triangulation_3<EPIC>  Cls;
  typedef CGAL::Delaunay_triangulation_3<EPEC>  Cls_with_epec;

  _test_cls_delaunay_3( Cls() );
  _test_cls_delaunay_3( Cls_with_epec() );

  typedef CGAL::Triangulation_data_structure_3<
    CGAL::Triangulation_vertex_base_3<K>,
    CGAL::Delaunay_triangulation_cell_base_3<K> >   Tds_Delaunay_Cb;
  typedef CGAL::Delaunay_triangulation_3<
    EPIC, Tds_Delaunay_Cb>                          Cls_with_Delaunay_Cb;

  _test_cls_delaunay_3( Cls_with_Delaunay_Cb() );

#ifdef CGAL_LINKED_WITH_TBB
  typedef CGAL::Spatial_lock_grid_3<
    CGAL::Tag_priority_blocking>                      Lock_ds;
  typedef CGAL::Triangulation_data_structure_3<
    CGAL::Triangulation_vertex_base_3<EPIC>,
    CGAL::Delaunay_triangulation_cell_base_3<EPIC>,
    CGAL::Parallel_tag >                                    Tds_parallel;
  typedef CGAL::Delaunay_triangulation_3<
    EPIC, Tds_parallel, CGAL::Default, Lock_ds>       Cls_parallel;
  // The following test won't do things in parallel since it doesn't provide
  // a lock data structure
  _test_cls_delaunay_3( Cls_parallel() );
  // This test performs parallel operations
  _test_cls_parallel_triangulation_3( Cls_parallel() );
#endif

  // Second version for the circumcenter storing cell base class.
  typedef CGAL::Triangulation_vertex_base_3<K>                          Vb;
  typedef CGAL::Delaunay_triangulation_cell_base_with_circumcenter_3<K> Cb;
  typedef CGAL::Triangulation_data_structure_3<Vb, Cb>                  TDS;
  typedef CGAL::Delaunay_triangulation_3<K, TDS>                        Cls_circumcenter;

  _test_cls_delaunay_3( Cls_circumcenter() );

  return 0;
}
